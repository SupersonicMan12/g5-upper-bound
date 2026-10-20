# Certificate export for an lubb_search.jl run (for the independent checker in independent/).
#
#   nice -n 19 julia --startup-file=no --project=. lubb_export.jl --logdir=4.5/logs_K4p5 --K=4.5 \
#         --output=4.5/certificates_K4p5.jsonl [--workers=6]
#
# For every logged node this script re-derives, exactly as lubb_replay.jl does, one
# certificate (nonnegative row multipliers) for
#   * every bound tightening (certified against the relaxation of the node's PRE-box), and
#   * every pruned leaf (certified on the stored, tightened box),
# checks each certificate in exact dyadic arithmetic (LUBB.certified_upper with Dy), and
# writes it out.  HiGHS only proposes multipliers; a certificate is written only after the
# exact check passes.  Rows are identified by NAME (see rowname below), not by position, so
# the file does not depend on the row order of lubb.jl.
#
# Output: JSON lines, one certificate per line.  Format documented in
# independent/code/data_formats.md ("Certificate export").
using Distributed, Dates, Printf, SHA

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
const KBITS = parse(Int, get(OPTS, "kbits", "20"))
const K = Float64(floor(KDEC * big(2)^KBITS) / big(2)^KBITS)
const LOGDIR = OPTS["logdir"]
const NWORK = parse(Int, get(OPTS, "workers", "6"))
const OUTPUT = OPTS["output"]
ispath(OUTPUT) && error("refusing to overwrite $OUTPUT")
isfile(joinpath(LOGDIR, "options.toml")) && error("this exporter handles version-1 logs only (no options.toml)")

NWORK > 0 && addprocs(NWORK; exeflags=["--project=$(@__DIR__)", "--startup-file=no", "--heap-size-hint=1G"])
@everywhere include(joinpath($(@__DIR__), "lubb.jl"))
@everywhere using .LUBB

@everywhere begin
    const B = LUBB
    const ST = B.build_structure(use_theorems=true, symmetry=false)
    const R = B.Relax(ST)
    const TOL_LADDER = [(1e-10, 4), (1e-10, 1), (1e-11, 4), (1e-11, 1), (1e-12, 1)]
    function settol!(tol)
        B.HiGHS.Highs_setDoubleOptionValue(R.h, "primal_feasibility_tolerance", tol)
        B.HiGHS.Highs_setDoubleOptionValue(R.h, "dual_feasibility_tolerance", tol)
    end

    # ---- row names, independent of the row order of lubb.jl ----------------------------
    const CUTRATIOS = (0.5, 0.75, 1.0, 1.25, 1.375, 1.5, 1.625, 1.75, 2.0)
    function rownames(st::B.Structure)
        names = String[]
        for k in 1:B.NN-1, i in k:B.NN, j in k+1:B.NN, sgn in ("+", "-")
            push!(names, "CP($k,$i,$j,$sgn)")          # sgn*A^(k)[i,j] <= p_k
        end
        for (i, t, j) in B.PIDX
            j == t || continue
            for sgn in ("+", "-")
                push!(names, "COL($i,$t,$sgn)")          # sgn*w[i,t,t] <= p_t
            end
        end
        for k in 1:B.NN-1
            push!(names, "DBL($k)")                       # p_{k+1} <= 2 p_k   (k=1: p_2 <= 2)
        end
        for k in 1:B.NN-2, r in CUTRATIOS
            push!(names, "COHEN($k,$r)")                  # p_{k+2} <= r^2 p_k + (3-2r) p_{k+1}
        end
        push!(names, "G4")                                # p_5 <= 4 p_2
        for (i, t, j) in B.PIDX
            for s in ("lolo", "hihi", "lohi", "hilo")
                push!(names, "MC($i,$t,$j,$s)")
            end
            if j > t
                for sgn in ("+", "-")
                    push!(names, "PERSP($i,$t,$j,$sgn)")
                end
            end
        end
        length(names) == length(st.rows) || error("row-name count $(length(names)) != rows $(length(st.rows))")
        return names
    end
    const ROWNAME = rownames(ST)
    # sanity: the k-th name must describe the k-th template
    function check_names()
        for (r, rt) in enumerate(ST.rows)
            nm = ROWNAME[r]
            if startswith(nm, "MC")
                rt.kind == 2 || error("name mismatch at row $r: $nm")
            elseif startswith(nm, "PERSP")
                rt.kind == 3 || error("name mismatch at row $r: $nm")
            else
                rt.kind == 1 || error("name mismatch at row $r: $nm")
            end
        end
    end
    check_names()

    # ---- certificates -------------------------------------------------------------------
    # A certificate: multipliers (sparse, by row name), objective (0 or ±e_v), and the exact
    # value Q of the Lemma 4.1 right side computed in Dy, plus what Q must satisfy.
    struct Cert
        lam::Vector{Tuple{Int,Float64}}   # (row index, multiplier > 0)
        objvar::Int                        # 0 for c = 0
        objcoef::Float64
        Q::Rational{BigInt}
    end
    sparse(lambda) = [(r, lambda[r]) for r in 1:length(lambda) if lambda[r] > 0]
    function try_lambda(box, y, cobj)
        for sgn in B.sign_order(y)
            lambda = [max(0.0, sgn * y[r]) for r in 1:R.nrow]
            Q = B.certified_upper(ST, box, lambda, cobj, B.Dy)
            Q === Inf && continue
            yield_ = (lambda, Rational{BigInt}(Q))
            return yield_
        end
        return nothing
    end
    # infeasibility certificate (c = 0, Q < 0) on `box` from the last solve / elastic LP
    function infeasible_cert(box)
        res = try_lambda(box, R.ray, Pair{Int,Float64}[])
        res !== nothing && res[2] < 0 && return Cert(sparse(res[1]), 0, 0.0, res[2])
        B.HiGHS.Highs_changeColsBoundsByRange(R.h, B.NVAR, B.NVAR, [0.0], [Inf])
        B.set_objective!(R, B.NVAR + 1, -1)
        s = B.solve!(R)
        out = nothing
        if s == :optimal
            res = try_lambda(box, R.rowdual, Pair{Int,Float64}[])
            res !== nothing && res[2] < 0 && (out = Cert(sparse(res[1]), 0, 0.0, res[2]))
        end
        B.HiGHS.Highs_changeColsBoundsByRange(R.h, B.NVAR, B.NVAR, [0.0], [0.0])
        return out
    end
    # one attempt at: sense*x_v <= target on R(box).  Returns Cert or nothing.
    function bound_attempt(box, v, sense, target)
        B.set_objective!(R, v, sense)
        s = B.solve!(R)
        if s == :infeasible
            c = infeasible_cert(box)
            c !== nothing && return c
        elseif s == :optimal
            res = try_lambda(box, R.rowdual, [v => Float64(sense)])
            res !== nothing && res[2] <= Rational{BigInt}(target) &&
                return Cert(sparse(res[1]), v, Float64(sense), res[2])
        end
        return nothing
    end
    function bound_cert(box, v, sense, target)
        B.set_box!(R, box); B.setstrat!(R, 4)
        c = bound_attempt(box, v, sense, target); c !== nothing && return c
        B.HiGHS.Highs_clearSolver(R.h)
        c = bound_attempt(box, v, sense, target); c !== nothing && return c
        for (tol, strat) in TOL_LADDER
            settol!(tol); B.setstrat!(R, strat); B.HiGHS.Highs_clearSolver(R.h)
            c = bound_attempt(box, v, sense, target)
            c !== nothing && break
        end
        settol!(1e-9); B.setstrat!(R, 4); B.HiGHS.Highs_clearSolver(R.h)
        return c
    end
    # leaf: no point of R(box) has p5 >= K.  Either c=0, Q<0, or c=e_p5, Q<K.
    function leaf_attempt(box)
        B.set_box!(R, box); B.set_objective!(R, B.P5, 1); B.setstrat!(R, 1)
        s = B.solve!(R)
        if s == :infeasible
            c = infeasible_cert(box); c !== nothing && return c
        elseif s == :optimal
            res = try_lambda(box, R.rowdual, [B.P5 => 1.0])
            res !== nothing && res[2] < Rational{BigInt}(K) && return Cert(sparse(res[1]), B.P5, 1.0, res[2])
        end
        # max p5 without its lower bound K; multipliers are evaluated on the actual box
        B.set_objective!(R, B.P5, 1)
        B.HiGHS.Highs_changeColsBoundsByRange(R.h, B.P5 - 1, B.P5 - 1, [-Inf], [box.hi[B.P5]])
        s = B.solve!(R)
        out = nothing
        if s == :optimal
            res = try_lambda(box, R.rowdual, [B.P5 => 1.0])
            res !== nothing && res[2] < Rational{BigInt}(K) && (out = Cert(sparse(res[1]), B.P5, 1.0, res[2]))
        end
        B.HiGHS.Highs_changeColsBoundsByRange(R.h, B.P5 - 1, B.P5 - 1, [box.lo[B.P5]], [box.hi[B.P5]])
        out !== nothing && return out
        B.set_box!(R, box); B.set_objective!(R, B.P5, 1); B.solve!(R)
        return infeasible_cert(box)
    end
    function leaf_cert(box)
        c = leaf_attempt(box); c !== nothing && return c
        for (tol, strat) in TOL_LADDER
            settol!(tol); B.HiGHS.Highs_clearSolver(R.h)
            c = leaf_attempt(box)
            c !== nothing && break
        end
        settol!(1e-9); B.HiGHS.Highs_clearSolver(R.h)
        return c
    end

    # ---- JSON output --------------------------------------------------------------------
    jnum(x::Float64) = repr(x)                 # shortest round-trip decimal: exact double
    function jcert(io, id, kind, var, side, target, c::Cert)
        print(io, "{\"node\":", id, ",\"kind\":\"", kind, "\"")
        if kind == "tighten"
            print(io, ",\"var\":\"", B.varname(var), "\",\"side\":\"", side, "\",\"target\":", jnum(target))
        end
        if c.objvar == 0
            print(io, ",\"objective\":null,\"check\":\"Q<0\"")
        else
            print(io, ",\"objective\":{\"var\":\"", B.varname(c.objvar), "\",\"coef\":", jnum(c.objcoef), "}")
            print(io, kind == "leaf" ? ",\"check\":\"Q<K\"" : ",\"check\":\"Q<=target\"")
        end
        print(io, ",\"Q\":\"", numerator(c.Q), "/", denominator(c.Q), "\",\"lam\":[")
        first = true
        for (r, l) in c.lam
            first || print(io, ",")
            print(io, "[\"", ROWNAME[r], "\",", jnum(l), "]")
            first = false
        end
        println(io, "]}")
    end

    function child_box(parent, side)
        bv = parent.branchvar
        (parent.box.lo[bv] < parent.split < parent.box.hi[bv]) || return nothing
        b = copy(parent.box)
        side == 1 ? (b.hi[bv] = parent.split) : (b.lo[bv] = parent.split)
        return b
    end
    # export all records of one file; returns (text, fails, ncert_tighten, ncert_leaf)
    function export_file(path, external::Dict{Int64,B.Box})
        recs = B.read_records(path)
        local_ = Dict(r.id => r for r in recs)
        io = IOBuffer()
        fails = Tuple{Int64,String}[]
        nt = 0; nl = 0
        for r in recs
            pre = haskey(external, r.id) ? external[r.id] :
                  haskey(local_, r.parent) ? child_box(local_[r.parent], r.side) : nothing
            pre === nothing && (push!(fails, (r.id, "no pre-box")); continue)
            post = r.box
            for v in B.OBBT_VARS, sense in (-1, 1)
                tighter = sense == 1 ? post.hi[v] < pre.hi[v] : post.lo[v] > pre.lo[v]
                tighter || continue
                target = sense == 1 ? post.hi[v] : -post.lo[v]
                c = bound_cert(pre, v, sense, target)
                if c === nothing
                    push!(fails, (r.id, "tightening $(B.varname(v)) $(sense == 1 ? "hi" : "lo")"))
                else
                    jcert(io, r.id, "tighten", v, sense == 1 ? "hi" : "lo", target, c); nt += 1
                end
            end
            if r.code == 0
                c = leaf_cert(post)
                if c === nothing
                    push!(fails, (r.id, "leaf"))
                else
                    jcert(io, r.id, "leaf", 0, "", 0.0, c); nl += 1
                end
            end
        end
        return String(take!(io)), fails, nt, nl
    end
end

function main()
    t0 = time()
    files = sort(filter(f -> endswith(f, ".bin"), readdir(LOGDIR; join=true)))
    allrecs = Dict{Int64,Any}()
    owner = Dict{Int64,Int}()
    for (fi, f) in enumerate(files), r in LUBB.read_records(f)
        haskey(allrecs, r.id) && error("duplicate id $(r.id)")
        allrecs[r.id] = r; owner[r.id] = fi
    end
    roots = [id for (id, r) in allrecs if r.parent == -1]
    length(roots) == 1 || error("expected one root")
    external = [Dict{Int64,LUBB.Box}() for _ in files]
    external[owner[roots[1]]][roots[1]] = LUBB.initial_box(K; use_theorems=true)
    for (id, r) in allrecs
        r.parent == -1 && continue
        p = allrecs[r.parent]
        owner[r.parent] == owner[id] && continue
        bv = p.branchvar
        (p.box.lo[bv] < p.split < p.box.hi[bv]) || error("bad split at $(p.id)")
        b = copy(p.box)
        r.side == 1 ? (b.hi[bv] = p.split) : (b.lo[bv] = p.split)
        external[owner[id]][id] = b
    end
    n = length(allrecs)
    @printf("loaded %d records from %d files (%.1fs)\n", n, length(files), time() - t0); flush(stdout)
    jobs = [(files[fi], external[fi]) for fi in eachindex(files)]
    results = pmap(j -> export_file(j[1], j[2]), jobs; batch_size=1)
    fails = reduce(vcat, [r[2] for r in results]; init=Tuple{Int64,String}[])
    nt = sum(r[3] for r in results); nl = sum(r[4] for r in results)
    open(OUTPUT, "w") do io
        println(io, "{\"format\":\"N5LUBBCertificates-v1\",\"K\":", repr(K), ",\"K_rational\":\"",
                numerator(Rational{BigInt}(K)), "/", denominator(Rational{BigInt}(K)),
                "\",\"logdir\":\"", LOGDIR, "\",\"nodes\":", n, ",\"tightenings\":", nt,
                ",\"leaves\":", nl, ",\"failures\":", length(fails),
                ",\"lubb_sha256\":\"", bytes2hex(sha256(read(joinpath(@__DIR__, "lubb.jl")))),
                "\",\"date_utc\":\"", now(UTC), "\"}")
        for r in results
            print(io, r[1])
        end
    end
    @printf("%s: %d tightening + %d leaf certificates, %d failures, %.1fs -> %s\n",
            isempty(fails) ? "EXPORT OK" : "EXPORT INCOMPLETE", nt, nl, length(fails), time() - t0, OUTPUT)
    for (id, msg) in fails[1:min(end, 20)]
        println("  node $id: $msg")
    end
    return isempty(fails)
end

main() || exit(1)
