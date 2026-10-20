# Parallel driver for the LU branch-and-bound (lubb.jl).
#
#   julia --startup-file=no --project=. lubb_search.jl --K=4.5 --workers=8 \
#         --output=lubb_K4p5.toml [--logdir=lubb_logs_K4p5] [--budget=2000] [--no-theorems]
#
# The target K is rounded DOWN to a multiple of 2^-20, so that a completed search proves
# p5 < K_dyadic <= K for every normalized completely pivoted 5x5 matrix.
# A report with status = "proved" requires: zero unresolved leaves, zero open nodes, and a
# complete tree (sum of 2^-depth over leaves equal to 1 in exact arithmetic).
using Distributed, Dates, Printf, SHA, TOML

function parse_args(args)
    opts = Dict{String,String}()
    for a in args
        startswith(a, "--") || error("bad argument $a")
        if occursin('=', a)
            k, v = split(a[3:end], '='; limit=2); opts[k] = v
        else
            opts[a[3:end]] = "true"
        end
    end
    return opts
end

const OPTS = parse_args(ARGS)
const KDEC = parse(BigFloat, OPTS["K"])
const K = Float64(floor(KDEC * big(2)^20) / big(2)^20)
const NWORK = parse(Int, get(OPTS, "workers", "8"))
const BUDGET = parse(Int, get(OPTS, "budget", "2000"))
const OUTPUT = get(OPTS, "output", "")
const LOGDIR = get(OPTS, "logdir", "")
const THEOREMS = !haskey(OPTS, "no-theorems")
const MAXSEC = parse(Float64, get(OPTS, "max-seconds", "1e9"))
const FRONTIER = parse(Int, get(OPTS, "frontier", string(8 * NWORK)))

isempty(OUTPUT) || !ispath(OUTPUT) || error("refusing to overwrite $OUTPUT")
isempty(LOGDIR) || mkpath(LOGDIR)
const SRC = ["lubb.jl", "lubb_search.jl"]
const HASHES = Dict(f => bytes2hex(sha256(read(joinpath(@__DIR__, f)))) for f in SRC)

addprocs(NWORK; exeflags=["--project=$(@__DIR__)", "--startup-file=no"])
@everywhere include(joinpath($(@__DIR__), "lubb.jl"))
@everywhere using .LUBB

@everywhere function lubb_worker(jobs, results, K, theorems, logdir)
    st = LUBB.build_structure(use_theorems=theorems)
    R = LUBB.Relax(st)
    root = LUBB.initial_box(K; use_theorems=theorems)
    while true
        job = take!(jobs)
        job === nothing && break
        task, node, budget = job
        log = isempty(logdir) ? nothing : open(joinpath(logdir, "task_$(task).bin"), "w")
        r0 = R.lpcount; f0 = R.lpfail
        res = LUBB.run_task!(R, node, task, K, root; budget, log)
        log === nothing || close(log)
        put!(results, (task, res.nodes, res.pruned, res.unresolved, res.maxdepth, res.open,
                       R.lpcount - r0, R.lpfail - f0, myid()))
    end
end

function main()
    t0 = time()
    st = LUBB.build_structure(use_theorems=THEOREMS)
    root = LUBB.initial_box(K; use_theorems=THEOREMS)
    # serial expansion of the first levels (task 0), recorded like any other task
    R = LUBB.Relax(st)
    log0 = isempty(LOGDIR) ? nothing : open(joinpath(LOGDIR, "task_0.bin"), "w")
    frontier = [LUBB.TaskNode(root, 0, Int64(-1), Int8(0))]
    counter = Int64(0)
    nodes = 0; unresolved = Tuple{Int64,LUBB.Box}[]; maxdepth = 0; lps = 0; lpfail = 0
    while !isempty(frontier) && length(frontier) < FRONTIER
        n = popfirst!(frontier)
        counter += 1
        id = counter
        res = LUBB.process_node!(R, n.box, K, root)
        nodes += 1
        log0 === nothing || LUBB.write_task_record(log0, id, n, res)
        if res.outcome == :branch
            q = res.branchvar
            b1 = copy(res.box); b1.hi[q] = res.split
            b2 = copy(res.box); b2.lo[q] = res.split
            push!(frontier, LUBB.TaskNode(b1, n.depth + 1, id, Int8(1)))
            push!(frontier, LUBB.TaskNode(b2, n.depth + 1, id, Int8(2)))
        elseif res.outcome == :unresolved
            push!(unresolved, (id, res.box))
        end
    end
    log0 === nothing || close(log0)
    lps += R.lpcount; lpfail += R.lpfail
    @printf("K=%.10f (dyadic %s) serial frontier=%d after %d nodes, %.1fs\n", K, string(K), length(frontier), nodes, time() - t0)
    flush(stdout)

    jobs = RemoteChannel(() -> Channel{Any}(10^7))
    results = RemoteChannel(() -> Channel{Any}(10^7))
    for p in workers()
        remote_do(lubb_worker, p, jobs, results, K, THEOREMS, LOGDIR)
    end
    nexttask = Int64(1)
    outstanding = 0
    for n in frontier
        put!(jobs, (nexttask, n, BUDGET)); nexttask += 1; outstanding += 1
    end
    lastreport = time()
    stopped = false
    open_left = LUBB.TaskNode[]
    while outstanding > 0
        task, tn, tp, tu, tmd, topen, tlps, tfail, wid = take!(results)
        outstanding -= 1
        nodes += tn; append!(unresolved, tu); maxdepth = max(maxdepth, tmd)
        lps += tlps; lpfail += tfail
        if time() - t0 > MAXSEC || length(unresolved) > 50
            stopped = true
        end
        for n in topen
            if stopped
                push!(open_left, n)
            else
                put!(jobs, (nexttask, n, BUDGET)); nexttask += 1; outstanding += 1
            end
        end
        if time() - lastreport > 60
            @printf("  %s nodes=%d outstanding=%d unresolved=%d maxdepth=%d lps=%d fail=%d %.0fs\n",
                    Dates.format(now(), "HH:MM:SS"), nodes, outstanding, length(unresolved), maxdepth,
                    lps, lpfail, time() - t0)
            flush(stdout)
            lastreport = time()
        end
    end
    for p in workers()
        put!(jobs, nothing)
    end
    elapsed = time() - t0
    # sources must not have changed during the run
    for f in SRC
        bytes2hex(sha256(read(joinpath(@__DIR__, f)))) == HASHES[f] || error("source $f changed during the run")
    end
    proved = isempty(unresolved) && isempty(open_left) && !stopped
    status = proved ? "proved" : "unresolved"
    @printf("%s: K=%s nodes=%d unresolved=%d open=%d maxdepth=%d lps=%d lpfail=%d elapsed=%.1fs\n",
            status, string(K), nodes, length(unresolved), length(open_left), maxdepth, lps, lpfail, elapsed)
    if !isempty(OUTPUT)
        rep = Dict{String,Any}(
            "format" => "N5LUBB-v1",
            "status" => status,
            "statement" => "If status = proved (and the replay confirms it), no real 5x5 matrix with any legal complete-pivoting sequence has fifth pivot p5 >= K; with g4 = 4 this gives g5 < K.",
            "K_requested" => OPTS["K"],
            "K_dyadic" => string(K),
            "K_rational" => string(Rational{BigInt}(K)),
            "use_theorems" => THEOREMS,
            "theorem_inputs" => THEOREMS ? ["g3 = 9/4 (Tornheim 1965, Cohen 1974)", "g4 = 4 (Cryer 1968, Cohen 1974)", "Cohen 1974: p3 <= 3 p2 - p2^2 for normalized 3x3 CP matrices with p3 > 2"] : ["none beyond the 2x2 bound p_{k+1} <= 2 p_k"],
            "nodes" => nodes, "lp_solves" => lps, "lp_failures" => lpfail,
            "max_depth" => maxdepth, "unresolved_count" => length(unresolved),
            "open_count" => length(open_left), "tasks" => nexttask - 1,
            "workers" => NWORK, "budget_per_task" => BUDGET, "elapsed_seconds" => elapsed,
            "date_utc" => string(now(UTC)), "julia_version" => string(VERSION),
            "code_sha256" => HASHES, "logdir" => LOGDIR,
            "arithmetic" => "HiGHS floating-point LPs suggest multipliers; every pruning and bound tightening is certified by outward-rounded interval arithmetic on exact dyadic box data",
            "unresolved_boxes" => [Dict("id" => id, "lo" => b.lo, "hi" => b.hi) for (id, b) in unresolved[1:min(end, 50)]],
        )
        open(OUTPUT, "w") do io
            TOML.print(io, rep; sorted=true)
        end
    end
    return proved
end

main()
