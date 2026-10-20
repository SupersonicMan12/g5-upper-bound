# LUBB: exhaustive branch-and-bound over the LU factors of a completely pivoted 5x5 matrix.
# Version 2 (2026-10-01): sparse certificate evaluation, optional transpose symmetry and
# tightening of L coordinates, solve filtering.  Version 1 is archived in lubb_runs/code_v1/.
# Version 3 (2026-10-01): optional exclusion boxes (outcome :excluded, log code 3), used
# together with the local lemma (local_lemma.jl) around the optimal matrix.
#
# Mathematical setting (see lubb_proof.md).  After permuting rows/columns into the complete
# pivoting order, scaling so that max|a_ij| = a_11 = 1, and applying the sign changes
# A -> D1*A*D2, every real 5x5 matrix with a legal complete-pivoting sequence (ties allowed)
# has an LU factorization A = L*U with
#   * L unit lower triangular, L[i,1] in [0,1] (i >= 2), L[i,t] in [-1,1] (t >= 2);
#   * U upper triangular, U[1,1] = 1, pivots U[t,t] = p_t > 0;
#   * for 1 <= k <= 4 and i,j >= k:  |sum_{t=k}^{min(i,j)} L[i,t]*U[t,j]| <= U[k,k]
#     (the (i,j) entry of the k-th Schur complement is bounded by the k-th pivot).
# The growth factor is max_k p_k.  To show that p_5 >= K is impossible, the search covers
# a box of (L, p, U_offdiagonal) values and proves every leaf box infeasible by an explicit
# nonnegative combination of valid linear inequalities (a Farkas certificate), checked in
# outward-rounded interval arithmetic (and, in the replay, in exact rational arithmetic).
#
# Products L[i,t]*U[t,j] are represented by auxiliary variables w with McCormick and
# "perspective" inequalities that are valid on the current box.  Floating-point LP solves
# (HiGHS) only suggest certificates; they are never trusted.

module LUBB

using HiGHS
using Printf

export Box, Iv, Relax, initial_box, process_node!, search_dfs!, NodeResult

# ---------------------------------------------------------------------------------------
# Outward-rounded interval arithmetic on Float64.  For +,-,* in round-to-nearest mode the
# exact result lies in [prevfloat(r), nextfloat(r)] where r is the computed value.
# ---------------------------------------------------------------------------------------
struct Iv
    lo::Float64
    hi::Float64
end
Iv(x::Float64) = Iv(x, x)
@inline dn(x::Float64) = prevfloat(x)
@inline up(x::Float64) = nextfloat(x)
Base.:+(a::Iv, b::Iv) = Iv(dn(a.lo + b.lo), up(a.hi + b.hi))
Base.:-(a::Iv, b::Iv) = Iv(dn(a.lo - b.hi), up(a.hi - b.lo))
Base.:-(a::Iv) = Iv(-a.hi, -a.lo)
function Base.:*(a::Iv, b::Iv)
    p1 = a.lo * b.lo; p2 = a.lo * b.hi; p3 = a.hi * b.lo; p4 = a.hi * b.hi
    return Iv(dn(min(p1, p2, p3, p4)), up(max(p1, p2, p3, p4)))
end
Base.zero(::Type{Iv}) = Iv(0.0, 0.0)
Base.zero(::Iv) = Iv(0.0, 0.0)
half(x::Iv) = Iv(dn(x.lo / 2), up(x.hi / 2))    # exact unless subnormal; widen anyway
half(x::Float64) = x / 2
half(x::Rational) = x / 2
lower(x::Iv) = x.lo
upper(x::Iv) = x.hi
lower(x::Real) = x
upper(x::Real) = x
isfinite_iv(x::Iv) = isfinite(x.lo) && isfinite(x.hi)

# Exact dyadic numbers m * 2^e (all data of the relaxation are dyadic: Float64 box bounds,
# their products, halves and sums, and Float64 multipliers).  Used by the exact replay.
struct Dy <: Real
    m::BigInt
    e::Int
end
function Dy(x::Float64)
    isfinite(x) || throw(DomainError(x, "non-finite value in exact arithmetic"))
    x == 0 && return Dy(big(0), 0)
    fr, ex = frexp(x)
    return Dy(BigInt(ldexp(fr, 53)), ex - 53)
end
function Base.:+(a::Dy, b::Dy)
    a.m == 0 && return b
    b.m == 0 && return a
    a.e <= b.e ? Dy(a.m + (b.m << (b.e - a.e)), a.e) : Dy((a.m << (a.e - b.e)) + b.m, b.e)
end
Base.:-(a::Dy) = Dy(-a.m, a.e)
Base.:-(a::Dy, b::Dy) = a + (-b)
Base.:*(a::Dy, b::Dy) = Dy(a.m * b.m, a.e + b.e)
Base.zero(::Type{Dy}) = Dy(big(0), 0)
Base.zero(::Dy) = Dy(big(0), 0)
Base.sign(a::Dy) = sign(a.m)
Base.:<(a::Dy, b::Dy) = sign(b - a) > 0
Base.:<=(a::Dy, b::Dy) = sign(b - a) >= 0
Base.:(==)(a::Dy, b::Dy) = sign(a - b) == 0
Base.isless(a::Dy, b::Dy) = a < b
Base.max(a::Dy, b::Dy) = a < b ? b : a
Base.min(a::Dy, b::Dy) = a < b ? a : b
half(a::Dy) = Dy(a.m, a.e - 1)
Base.Rational{BigInt}(a::Dy) = a.e >= 0 ? Rational{BigInt}(a.m << a.e) : a.m // (big(1) << (-a.e))
Base.Float64(a::Dy) = Float64(Rational{BigInt}(a))

conv(::Type{Iv}, x::Float64) = Iv(x, x)
conv(::Type{Float64}, x::Float64) = x
conv(::Type{Rational{BigInt}}, x::Float64) = Rational{BigInt}(x)
conv(::Type{Dy}, x::Float64) = Dy(x)
conv(::Type{T}, x::Integer) where {T} = conv(T, Float64(x))

# ---------------------------------------------------------------------------------------
# Variables.  1:10 L (column-major strictly lower), 11:14 pivots p2..p5, 15:24 U strictly
# upper (row-major), 25:60 products w[i,t,j] = L[i,t]*U[t,j].
# ---------------------------------------------------------------------------------------
const NN = 5
const LIDX = [(i, t) for t in 1:NN-1 for i in t+1:NN]
const UOFF = [(t, j) for t in 1:NN-1 for j in t+1:NN]
const PIDX = [(i, t, j) for t in 1:NN-1 for i in t+1:NN for j in t:NN if !(t == 1 && j == 1)]
const NL = length(LIDX)
const NBOX = NL + 4 + length(UOFF)        # 24 branching coordinates
const NVAR = NBOX + length(PIDX)          # 60 LP variables
const LPOS = Dict(LIDX[q] => q for q in eachindex(LIDX))
const UPOS = Dict(UOFF[q] => NL + 4 + q for q in eachindex(UOFF))
const WPOS = Dict(PIDX[q] => NBOX + q for q in eachindex(PIDX))
Lvar(i, t) = LPOS[(i, t)]
function Uvar(t, j)
    t == 1 && j == 1 && return 0              # the constant U[1,1] = 1
    t == j && return NL + t - 1               # pivot p_t
    return UPOS[(t, j)]
end
Wvar(i, t, j) = WPOS[(i, t, j)]
const PIVOT_VAR = [0, NL + 1, NL + 2, NL + 3, NL + 4]
const P5 = NL + 4

function varname(v)
    v <= NL && return "L$(LIDX[v][1])$(LIDX[v][2])"
    v <= NL + 4 && return "p$(v - NL + 1)"
    v <= NBOX && return "U$(UOFF[v-NL-4][1])$(UOFF[v-NL-4][2])"
    (i, t, j) = PIDX[v-NBOX]
    return "w$(i)$(t)$(j)"
end

# ---------------------------------------------------------------------------------------
# Row templates.  Every row reads  sum_k coef_k * x[var_k] <= rhs.
# kind 1: constant row; kind 2: McCormick (product q, case s); kind 3: perspective (q, sign)
# ---------------------------------------------------------------------------------------
struct RowTemplate
    vars::Vector{Int}
    coefs::Vector{Float64}   # used for constant rows only
    rhs::Float64
    kind::Int
    q::Int
    s::Int
end

struct Structure
    rows::Vector{RowTemplate}
    use_theorems::Bool
    symmetry::Bool
end

"""
    build_structure(; use_theorems=true, symmetry=false, cutratios=...)

Constant rows:
  * complete-pivoting rows  ±A^(k)[i,j] - p_k <= 0  (k = 1: the pivot is the constant 1);
  * column rows  ±w[i,t,t] - p_t <= 0   (|L[i,t]| <= 1 times p_t);
  * elementary pivot rows  p_{k+1} <= 2 p_k (2x2 growth);
  * if use_theorems: Cohen's three-pivot inequality in tangent form
      p_{k+2} <= r^2 p_k + (3 - 2r) p_{k+1}   for dyadic r,
    and p_5 <= 4 p_2 (Cryer's g_4 = 4 applied to the second Schur complement);
  * if symmetry: |U12| <= L21.  Transposition maps a normalized CP matrix to one whose
    (L21, |U12|) are swapped (same pivots), so this is no loss of generality.
"""
function build_structure(; use_theorems::Bool=true, symmetry::Bool=false,
                         cutratios=(0.5, 0.75, 1.0, 1.25, 1.375, 1.5, 1.625, 1.75, 2.0))
    rows = RowTemplate[]
    if symmetry
        push!(rows, RowTemplate([Uvar(1, 2), Lvar(2, 1)], [1.0, -1.0], 0.0, 1, 0, 0))
        push!(rows, RowTemplate([Uvar(1, 2), Lvar(2, 1)], [-1.0, -1.0], 0.0, 1, 0, 0))
    end
    term(i, t, j) = t == i ? Uvar(i, j) : Wvar(i, t, j)
    for k in 1:NN-1, i in k:NN, j in k+1:NN
        vs = [term(i, t, j) for t in k:min(i, j)]
        for sgn in (1.0, -1.0)
            vars = copy(vs); coefs = fill(sgn, length(vs)); rhs = 0.0
            if k == 1
                rhs = 1.0
            else
                push!(vars, PIVOT_VAR[k]); push!(coefs, -1.0)
            end
            push!(rows, RowTemplate(vars, coefs, rhs, 1, 0, 0))
        end
    end
    for (q, (i, t, j)) in enumerate(PIDX)
        if j == t
            for sgn in (1.0, -1.0)
                push!(rows, RowTemplate([Wvar(i, t, j), PIVOT_VAR[t]], [sgn, -1.0], 0.0, 1, 0, 0))
            end
        end
    end
    for k in 1:NN-1
        if k == 1
            push!(rows, RowTemplate([PIVOT_VAR[2]], [1.0], 2.0, 1, 0, 0))
        else
            push!(rows, RowTemplate([PIVOT_VAR[k+1], PIVOT_VAR[k]], [1.0, -2.0], 0.0, 1, 0, 0))
        end
    end
    if use_theorems
        for k in 1:NN-2, r in cutratios
            @assert r * 64 == round(r * 64)   # dyadic, so r^2 and 3-2r are exact
            if k == 1
                push!(rows, RowTemplate([PIVOT_VAR[3], PIVOT_VAR[2]], [1.0, -(3 - 2r)], r^2, 1, 0, 0))
            else
                push!(rows, RowTemplate([PIVOT_VAR[k+2], PIVOT_VAR[k], PIVOT_VAR[k+1]],
                                        [1.0, -r^2, -(3 - 2r)], 0.0, 1, 0, 0))
            end
        end
        push!(rows, RowTemplate([P5, PIVOT_VAR[2]], [1.0, -4.0], 0.0, 1, 0, 0))
    end
    # drop explicit zero coefficients from constant rows (e.g. r = 3/2 gives 3 - 2r = 0)
    rows = [rt.kind == 1 ? RowTemplate(rt.vars[rt.coefs .!= 0], rt.coefs[rt.coefs .!= 0],
                                       rt.rhs, 1, 0, 0) : rt for rt in rows]
    for (q, (i, t, j)) in enumerate(PIDX)
        w = Wvar(i, t, j); L = Lvar(i, t); U = Uvar(t, j)
        for s in 1:4
            push!(rows, RowTemplate([w, L, U], Float64[], 0.0, 2, q, s))
        end
        if j > t
            for s in 1:2
                vars = t == 1 ? [w, U] : [w, U, PIVOT_VAR[t]]
                push!(rows, RowTemplate(vars, Float64[], 0.0, 3, q, s))
            end
        end
    end
    return Structure(rows, use_theorems, symmetry)
end

# Box: lower/upper bounds of the 24 branching coordinates (Float64 = exact dyadic numbers).
struct Box
    lo::Vector{Float64}
    hi::Vector{Float64}
end
Base.copy(b::Box) = Box(copy(b.lo), copy(b.hi))

"""
    initial_box(K; use_theorems=true)

Root box under the hypothesis p_5 >= K (K >= 4).
Elementary facts: |L| <= 1 (L[:,1] >= 0 by normalization), |U[t,j]| <= p_t, p_2 <= 2,
p_{k+1} <= 2 p_k.  With use_theorems: g_3 = 9/4 and g_4 = 4 give p_3 <= 9/4, p_4 <= 4,
p_5 <= (9/4)^2, p_4 >= K/2, p_3 >= K/(9/4), p_2 >= K/4.  Without theorems only the
2x2 bound is used: p_3 <= 4, p_4 <= 8, p_5 <= 16, p_4 >= K/2, p_3 >= K/4, p_2 >= K/8.
Lower bounds are rounded down and upper bounds up, to dyadic numbers.
"""
function initial_box(K::Float64; use_theorems::Bool=true)
    lo = zeros(NBOX); hi = zeros(NBOX)
    for (q, (i, t)) in enumerate(LIDX)
        lo[q] = t == 1 ? 0.0 : -1.0; hi[q] = 1.0
    end
    if use_theorems
        plo = [K / 4, K / 2.25, K / 2, K]
        phi = [2.0, 2.25, 4.0, 5.0625]
    else
        plo = [K / 8, K / 4, K / 2, K]
        phi = [2.0, 4.0, 8.0, 16.0]
    end
    for t in 1:4
        lo[NL+t] = rounddown_dyadic(plo[t]); hi[NL+t] = phi[t]
    end
    lo[P5] = K   # the hypothesis p5 >= K itself (K must be a Float64, i.e. dyadic)
    for (q, (t, j)) in enumerate(UOFF)
        h = t == 1 ? 1.0 : hi[PIVOT_VAR[t]]
        lo[NL+4+q] = -h; hi[NL+4+q] = h
    end
    return Box(lo, hi)
end

# round x (a Float64 result of a division) down safely to a 2^-30 dyadic grid
rounddown_dyadic(x) = floor(prevfloat(x) * 2.0^30) / 2.0^30

# ---------------------------------------------------------------------------------------
# Row values for a box, generic over arithmetic T in {Float64, Iv, Rational{BigInt}}.
# Returns coefficient vectors aligned with template.vars and the right-hand sides.
# ---------------------------------------------------------------------------------------
function bnd(box::Box, v::Int, ::Type{T}) where {T}
    v == 0 && return (conv(T, 1.0), conv(T, 1.0))
    return (conv(T, box.lo[v]), conv(T, box.hi[v]))
end

# Coefficients of one row for a box (written into cr, aligned with rt.vars); returns the
# right-hand side.
function row_value!(cr::AbstractVector{T}, rt::RowTemplate, box::Box) where {T}
    if rt.kind == 1
        for k in eachindex(rt.vars)
            cr[k] = conv(T, rt.coefs[k])
        end
        return conv(T, rt.rhs)
    elseif rt.kind == 2
        (i, t, j) = PIDX[rt.q]
        a, b = bnd(box, Lvar(i, t), T)
        c, d = bnd(box, Uvar(t, j), T)
        one_ = conv(T, 1.0)
        if rt.s == 1      # -w + c L + a U <= a c
            cr[1] = -one_; cr[2] = c; cr[3] = a; return a * c
        elseif rt.s == 2  # -w + d L + b U <= b d
            cr[1] = -one_; cr[2] = d; cr[3] = b; return b * d
        elseif rt.s == 3  #  w - d L - a U <= -a d
            cr[1] = one_; cr[2] = -d; cr[3] = -a; return -(a * d)
        else              #  w - c L - b U <= -b c
            cr[1] = one_; cr[2] = -c; cr[3] = -b; return -(b * c)
        end
    else
        (i, t, j) = PIDX[rt.q]
        a, b = bnd(box, Lvar(i, t), T)
        m = half(a + b); rad = half(b - a)
        sg = rt.s == 1 ? 1.0 : -1.0
        cr[1] = conv(T, sg)                     # ±w
        cr[2] = sg == 1.0 ? -m : m              # ∓m U
        if t == 1                               # |U_1j| <= 1
            return rad
        else                                    # |U_tj| <= p_t
            cr[3] = -rad
            return conv(T, 0.0)
        end
    end
end

function row_values!(coefs::Vector{Vector{T}}, rhs::Vector{T}, st::Structure, box::Box) where {T}
    for (r, rt) in enumerate(st.rows)
        rhs[r] = row_value!(coefs[r], rt, box)
    end
    return coefs, rhs
end

function alloc_rows(st::Structure, ::Type{T}) where {T}
    coefs = [Vector{T}(undef, length(rt.vars)) for rt in st.rows]
    rhs = Vector{T}(undef, length(st.rows))
    return coefs, rhs
end

# bounds of every LP variable (box coordinates, then products from corner values)
function var_bounds(box::Box, ::Type{T}) where {T}
    lo = Vector{T}(undef, NVAR); hi = Vector{T}(undef, NVAR)
    for v in 1:NBOX
        lo[v] = conv(T, box.lo[v]); hi[v] = conv(T, box.hi[v])
    end
    for (q, (i, t, j)) in enumerate(PIDX)
        a, b = bnd(box, Lvar(i, t), T)
        c, d = bnd(box, Uvar(t, j), T)
        ps = (a * c, a * d, b * c, b * d)
        lo[NBOX+q] = conv_min(ps); hi[NBOX+q] = conv_max(ps)
    end
    return lo, hi
end
conv_min(ps::NTuple{4,Iv}) = Iv(minimum(p.lo for p in ps), minimum(p.lo for p in ps))
conv_max(ps::NTuple{4,Iv}) = Iv(maximum(p.hi for p in ps), maximum(p.hi for p in ps))
conv_min(ps::NTuple{4,T}) where {T<:Real} = minimum(ps)
conv_max(ps::NTuple{4,T}) where {T<:Real} = maximum(ps)

# ---------------------------------------------------------------------------------------
# Certificate evaluation.  For multipliers lambda >= 0 and objective c (sparse Dict), every
# point x of the box satisfying all rows has
#     c.x <= lambda.rhs + max_{x in box} (c - A^T lambda).x .
# With c = 0 a negative value proves infeasibility.  Returned values are rigorous upper
# bounds (interval upper endpoints, or exact rationals).
# ---------------------------------------------------------------------------------------
function certified_upper(st::Structure, box::Box, lambda::Vector{Float64},
                         cobj::Vector{Pair{Int,Float64}}, ::Type{T}) where {T}
    vlo, vhi = var_bounds(box, T)
    g = [conv(T, 0.0) for _ in 1:NVAR]
    beta = conv(T, 0.0)
    cr = Vector{T}(undef, 4)
    for (r, rt) in enumerate(st.rows)
        l = lambda[r]
        l > 0 || continue                # only the rows used by the certificate are formed
        isfinite(l) || return Inf
        length(cr) < length(rt.vars) && resize!(cr, length(rt.vars))
        rh = row_value!(cr, rt, box)
        L = conv(T, l)
        for k in eachindex(rt.vars)
            v = rt.vars[k]
            g[v] = g[v] + L * cr[k]
        end
        beta = beta + L * rh
    end
    for (v, cv) in cobj
        g[v] = g[v] - conv(T, cv)       # g now holds A^T lambda - c
    end
    total = beta
    for v in 1:NVAR
        d = -g[v]                        # d = c - A^T lambda
        # max over x in [vlo, vhi] of d*x
        total = total + maxprod(d, vlo[v], vhi[v])
    end
    return upper(total)
end

function maxprod(d::Iv, xlo::Iv, xhi::Iv)
    p1 = d * Iv(xlo.lo, xlo.lo); p2 = d * Iv(xhi.hi, xhi.hi)
    return Iv(max(p1.lo, p2.lo), max(p1.hi, p2.hi))
end
maxprod(d::T, xlo::T, xhi::T) where {T<:Real} = max(d * xlo, d * xhi)

# ---------------------------------------------------------------------------------------
# HiGHS model (floating point, untrusted oracle).
# ---------------------------------------------------------------------------------------
mutable struct Relax
    st::Structure
    h::Ptr{Cvoid}
    nrow::Int
    coefsF::Vector{Vector{Float64}}
    rhsF::Vector{Float64}
    colval::Vector{Float64}
    coldual::Vector{Float64}
    rowval::Vector{Float64}
    rowdual::Vector{Float64}
    ray::Vector{Float64}
    varrows::Vector{Int}      # rows with box-dependent coefficients
    lpcount::Int
    lpfail::Int
end

function Relax(st::Structure)
    h = HiGHS.Highs_create()
    HiGHS.Highs_setBoolOptionValue(h, "output_flag", false)
    HiGHS.Highs_setStringOptionValue(h, "presolve", "off")
    HiGHS.Highs_setStringOptionValue(h, "solver", "simplex")
    HiGHS.Highs_setDoubleOptionValue(h, "primal_feasibility_tolerance", 1e-9)
    HiGHS.Highs_setDoubleOptionValue(h, "dual_feasibility_tolerance", 1e-9)
    nrow = length(st.rows)
    coefsF, rhsF = alloc_rows(st, Float64)
    box0 = initial_box(4.5; use_theorems=st.use_theorems)
    row_values!(coefsF, rhsF, st, box0)
    # row-wise matrix; replace structural zeros by placeholders so the pattern is kept.
    # Column NVAR+1 is an elastic slack t (every row reads a.x - t <= b); it is fixed to 0
    # except when searching for a well-conditioned infeasibility certificate.
    astart = Cint[]; aindex = Cint[]; avalue = Float64[]
    for (r, rt) in enumerate(st.rows)
        push!(astart, length(aindex))
        for k in eachindex(rt.vars)
            push!(aindex, rt.vars[k] - 1)
            val = coefsF[r][k]
            # placeholders only for box-dependent entries (overwritten by set_box!)
            push!(avalue, (rt.kind != 1 && val == 0.0) ? 1.0 : val)
        end
        push!(aindex, NVAR); push!(avalue, -1.0)
    end
    cost = zeros(NVAR + 1); cost[P5] = 1.0
    vlo, vhi = var_bounds(box0, Float64)
    push!(vlo, 0.0); push!(vhi, 0.0)
    rlo = fill(-Inf, nrow)
    ret = HiGHS.Highs_passLp(h, NVAR + 1, nrow, length(aindex), HiGHS.kHighsMatrixFormatRowwise,
        HiGHS.kHighsObjSenseMaximize, 0.0, cost, vlo, vhi, rlo, rhsF, astart, aindex, avalue)
    ret == HiGHS.kHighsStatusError && error("Highs_passLp failed")
    varrows = [r for (r, rt) in enumerate(st.rows) if rt.kind != 1]
    return Relax(st, h, nrow, coefsF, rhsF, zeros(NVAR + 1), zeros(NVAR + 1), zeros(nrow),
                 zeros(nrow), zeros(nrow), varrows, 0, 0)
end

function set_box!(R::Relax, box::Box)
    row_values!(R.coefsF, R.rhsF, R.st, box)
    for r in R.varrows
        rt = R.st.rows[r]
        for k in eachindex(rt.vars)
            HiGHS.Highs_changeCoeff(R.h, r - 1, rt.vars[k] - 1, R.coefsF[r][k])
        end
    end
    HiGHS.Highs_changeRowsBoundsByRange(R.h, 0, R.nrow - 1, fill(-Inf, R.nrow), R.rhsF)
    vlo, vhi = var_bounds(box, Float64)
    HiGHS.Highs_changeColsBoundsByRange(R.h, 0, NVAR - 1, vlo, vhi)
    return nothing
end

function set_objective!(R::Relax, v::Int, sense::Int)
    cost = zeros(NVAR + 1); cost[v] = 1.0
    HiGHS.Highs_changeColsCostByRange(R.h, 0, NVAR, cost)
    HiGHS.Highs_changeObjectiveSense(R.h, sense == 1 ? HiGHS.kHighsObjSenseMaximize : HiGHS.kHighsObjSenseMinimize)
end

# Minimize the uniform row violation t.  A positive optimum proves infeasibility, and its
# row duals form a Farkas certificate whose margin is about t (well conditioned).
function elastic_infeasible!(R::Relax, box::Box)
    HiGHS.Highs_changeColsBoundsByRange(R.h, NVAR, NVAR, [0.0], [Inf])
    set_objective!(R, NVAR + 1, -1)
    s = solve!(R)
    ok = false
    if s == :optimal && R.colval[NVAR+1] > 0
        for sgn in sign_order(R.rowdual)
            lambda = [max(0.0, sgn * R.rowdual[r]) for r in 1:R.nrow]
            if certified_upper(R.st, box, lambda, Pair{Int,Float64}[], Iv) < 0
                ok = true; break
            end
        end
    end
    HiGHS.Highs_changeColsBoundsByRange(R.h, NVAR, NVAR, [0.0], [0.0])
    return ok
end

# returns :optimal, :infeasible or :failed
function solve!(R::Relax)
    for attempt in 1:2
        R.lpcount += 1
        HiGHS.Highs_run(R.h)
        st = HiGHS.Highs_getModelStatus(R.h)
        if st == HiGHS.kHighsModelStatusOptimal
            HiGHS.Highs_getSolution(R.h, R.colval, R.coldual, R.rowval, R.rowdual)
            return :optimal
        elseif st == HiGHS.kHighsModelStatusInfeasible
            has = Ref{Cint}(0)
            HiGHS.Highs_getDualRay(R.h, has, R.ray)
            has[] != 0 && return :infeasible
        end
        R.lpfail += 1
        HiGHS.Highs_clearSolver(R.h)
    end
    return :failed
end

# HiGHS simplex strategy: 1 = dual (after bound changes), 4 = primal (after objective changes)
setstrat!(R::Relax, s::Int) = HiGHS.Highs_setIntOptionValue(R.h, "simplex_strategy", s)

# sign order to try for a multiplier vector: the dominant sign first
function sign_order(y)
    pos = 0.0; neg = 0.0
    for x in y
        x > 0 ? (pos += x) : (neg -= x)
    end
    return neg >= pos ? (-1.0, 1.0) : (1.0, -1.0)
end

# Rigorous upper bound of `sense * x[v]` over the relaxation of `box`, from the last solve.
# Tries the multiplier sign conventions; any lambda >= 0 is mathematically valid.
function bound_from_last(R::Relax, box::Box, v::Int, sense::Int, status::Symbol, ::Type{T}=Iv) where {T}
    best = Inf
    y = status == :optimal ? R.rowdual : R.ray
    for sgn in sign_order(y)
        lambda = [max(0.0, sgn * y[r]) for r in 1:R.nrow]
        cobj = status == :optimal ? [v => Float64(sense)] : Pair{Int,Float64}[]
        ub = certified_upper(R.st, box, lambda, cobj, T)
        if status == :infeasible
            ub < 0 && return -Inf
        else
            best = min(best, ub)
            best < 1e6 && break          # the dominant sign gave a meaningful bound
        end
    end
    return status == :infeasible ? Inf : best
end

# ---------------------------------------------------------------------------------------
# Node processing
# ---------------------------------------------------------------------------------------
struct NodeResult
    outcome::Symbol          # :pruned, :branch, :unresolved
    box::Box                 # box after bound tightening
    branchvar::Int
    split::Float64
    lpvalue::Float64
end

const OBBT_VARS = vcat(collect(NL+1:NL+3), collect(NL+5:NBOX))   # p2..p4 and U off-diagonal
const OBBT_VARS_L = vcat(collect(1:NL), OBBT_VARS)                # ... and the L coordinates

"""
    process_node!(R, box, K; obbt=true)

Proves (rigorously) the box infeasible, or returns a tightened box and a branching choice.
Tightened bounds are only ever taken from certified values, relaxed outward by a tiny
margin so that an independent replay can re-certify them.
"""
# Certify infeasibility of the relaxation of `box` after a solve reported it: first the
# solver's Farkas ray, then (if that is too ill-conditioned) an elastic re-solve.
function certify_infeasible!(R::Relax, box::Box, v::Int, sense::Int)
    bound_from_last(R, box, v, sense, :infeasible) == -Inf && return true
    ok = elastic_infeasible!(R, box)
    set_objective!(R, v, sense)
    return ok
end

contained(b::Box, E::Box) = all(E.lo[v] <= b.lo[v] && b.hi[v] <= E.hi[v] for v in 1:NBOX)
overlaps(b::Box, E::Box) = all(b.lo[v] < E.hi[v] && E.lo[v] < b.hi[v] for v in 1:NBOX)

# split choice that separates the box from an exclusion box E it overlaps but is not inside
function carve_split(b::Box, E::Box)
    best = 0; bestext = -1.0; pt = 0.0
    for v in 1:NBOX
        w = b.hi[v] - b.lo[v]
        w > 0 || continue
        if b.lo[v] < E.lo[v] < b.hi[v] && (E.lo[v] - b.lo[v]) / w > bestext
            best = v; bestext = (E.lo[v] - b.lo[v]) / w; pt = E.lo[v]
        end
        if b.lo[v] < E.hi[v] < b.hi[v] && (b.hi[v] - E.hi[v]) / w > bestext
            best = v; bestext = (b.hi[v] - E.hi[v]) / w; pt = E.hi[v]
        end
    end
    return best, pt
end

function process_node!(R::Relax, box0::Box, K::Float64, root::Box; obbt::Bool=true,
                       obbt_vars=OBBT_VARS, minrel=1e-7, margin_rel=1e-5, margin_abs=1e-7,
                       filter::Bool=false, filter_tol=1e-3, excl::Vector{Box}=Box[])
    for (ie, E) in enumerate(excl)
        contained(box0, E) && return NodeResult(:excluded, box0, ie, 0.0, NaN)
    end
    box = copy(box0)
    set_box!(R, box)
    set_objective!(R, P5, 1)
    setstrat!(R, 1)
    st = solve!(R)
    if st == :infeasible
        certify_infeasible!(R, box, P5, 1) && return NodeResult(:pruned, box, 0, 0.0, -Inf)
    elseif st == :optimal
        bound_from_last(R, box, P5, 1, st) < K && return NodeResult(:pruned, box, 0, 0.0, -Inf)
    end
    if obbt && st != :failed
        # Every tightening is certified against the relaxation of the ORIGINAL node box
        # box0 (the LP that HiGHS solves), so the certificates do not depend on the order.
        # The certified bound is relaxed outward by margin_abs + margin_rel*width, leaving
        # room for an independent replay whose multipliers may differ slightly.
        # Filtering (a speed heuristic only): once some LP solution of this node has x_v
        # within filter_tol*width of hi_v, maximizing x_v could tighten hi_v by at most
        # that much, so the solve is skipped (likewise for lo_v).
        setstrat!(R, 4)
        skiphi = falses(NBOX); skiplo = falses(NBOX)
        function mark!(x)
            for u in 1:NBOX
                wu = box0.hi[u] - box0.lo[u]
                x[u] >= box0.hi[u] - filter_tol * wu && (skiphi[u] = true)
                x[u] <= box0.lo[u] + filter_tol * wu && (skiplo[u] = true)
            end
        end
        filter && st == :optimal && mark!(R.colval)
        for v in obbt_vars
            for sense in (-1, 1)
                filter && (sense == 1 ? skiphi[v] : skiplo[v]) && continue
                set_objective!(R, v, sense)
                s2 = solve!(R)
                filter && s2 == :optimal && mark!(R.colval)
                if s2 == :infeasible
                    if certify_infeasible!(R, box0, v, sense)
                        set_objective!(R, P5, 1)
                        return NodeResult(:pruned, box0, 0, 0.0, -Inf)
                    end
                elseif s2 == :optimal
                    ub = bound_from_last(R, box0, v, sense, s2)   # sense*x_v <= ub on box0
                    isfinite(ub) || continue
                    w = box0.hi[v] - box0.lo[v]
                    slack = margin_abs + margin_rel * w
                    if sense == 1
                        nb = ub + slack
                        nb < box.hi[v] && (box.hi[v] = max(nb, box.lo[v]))
                    else
                        nb = -ub - slack
                        nb > box.lo[v] && (box.lo[v] = min(nb, box.hi[v]))
                    end
                end
            end
        end
        set_objective!(R, P5, 1)
        set_box!(R, box)
        setstrat!(R, 1)
        st = solve!(R)
        if st == :infeasible
            certify_infeasible!(R, box, P5, 1) && return NodeResult(:pruned, box, 0, 0.0, -Inf)
        elseif st == :optimal
            bound_from_last(R, box, P5, 1, st) < K && return NodeResult(:pruned, box, 0, 0.0, -Inf)
        end
    end
    lpv = st == :optimal ? R.colval[P5] : NaN
    for (ie, E) in enumerate(excl)
        contained(box, E) && return NodeResult(:excluded, box, ie, 0.0, lpv)
    end
    for E in excl
        if overlaps(box, E)
            q, pt = carve_split(box, E)
            q != 0 && return NodeResult(:branch, box, q, pt, lpv)
        end
    end
    q = choose_branch(R, box, root, st == :optimal; minrel)
    q == 0 && return NodeResult(:unresolved, box, 0, 0.0, lpv)
    mid = (box.lo[q] + box.hi[q]) / 2
    (mid <= box.lo[q] || mid >= box.hi[q]) && return NodeResult(:unresolved, box, 0, 0.0, lpv)
    return NodeResult(:branch, box, q, mid, lpv)
end

function choose_branch(R::Relax, box::Box, root::Box, valid::Bool; minrel=1e-7)
    rel(q) = (box.hi[q] - box.lo[q]) / max(root.hi[q] - root.lo[q], 1e-300)
    score = zeros(NBOX)
    if valid
        x = R.colval
        for (qq, (i, t, j)) in enumerate(PIDX)
            ql = Lvar(i, t); qu = Uvar(t, j)
            viol = abs(x[NBOX+qq] - x[ql] * x[qu])
            (isfinite(viol) && viol < 1e3) || continue
            wl = box.hi[ql] - box.lo[ql]; wu = box.hi[qu] - box.lo[qu]
            ml = max(abs(box.lo[ql]), abs(box.hi[ql])); mu = max(abs(box.lo[qu]), abs(box.hi[qu]))
            tot = wl * mu + wu * ml + 1e-300
            score[ql] += viol * wl * mu / tot
            score[qu] += viol * wu * ml / tot
        end
    end
    best = 0; bs = -1.0
    for q in 1:NBOX
        q == P5 && continue
        rel(q) < minrel && continue
        s = score[q] + 1e-9 * rel(q)
        if s > bs
            bs = s; best = q
        end
    end
    if best != 0 && score[best] <= 1e-8
        # no informative violation: bisect the relatively widest coordinate
        best = 0; bs = -1.0
        for q in 1:NBOX
            q == P5 && continue
            if rel(q) >= minrel && rel(q) > bs
                bs = rel(q); best = q
            end
        end
    end
    return best
end

# ---------------------------------------------------------------------------------------
# Task-based depth-first search with a global node-id scheme for logging and replay.
# id = task * 2^32 + local counter.  side: 0 root, 1 lower half, 2 upper half of parent.
# Log record (little-endian): id Int64, parent Int64, side Int8, depth Int32, outcome Int8
# (0 pruned, 1 branch, 2 unresolved), branchvar Int32, split Float64, then the box after
# bound tightening as 24 (lo, hi) Float64 pairs.
# ---------------------------------------------------------------------------------------
struct TaskNode
    box::Box
    depth::Int
    parent::Int64
    side::Int8
end

const RECORD_BYTES = 8 + 8 + 1 + 4 + 1 + 4 + 8 + 16 * NBOX

function write_task_record(io, id::Int64, n::TaskNode, res::NodeResult)
    write(io, id, n.parent, n.side, Int32(n.depth))
    code = res.outcome == :pruned ? Int8(0) : res.outcome == :branch ? Int8(1) :
           res.outcome == :excluded ? Int8(3) : Int8(2)
    write(io, code, Int32(res.branchvar), res.split)
    for v in 1:NBOX
        write(io, res.box.lo[v], res.box.hi[v])
    end
end

function read_records(path)
    recs = NamedTuple[]
    open(path, "r") do io
        while !eof(io)
            id = read(io, Int64); parent = read(io, Int64); side = read(io, Int8)
            depth = read(io, Int32); code = read(io, Int8); bv = read(io, Int32); split = read(io, Float64)
            lo = zeros(NBOX); hi = zeros(NBOX)
            for v in 1:NBOX
                lo[v] = read(io, Float64); hi[v] = read(io, Float64)
            end
            push!(recs, (id=id, parent=parent, side=side, depth=Int(depth), code=code,
                         branchvar=Int(bv), split=split, box=Box(lo, hi)))
        end
    end
    return recs
end

function run_task!(R::Relax, start::TaskNode, task::Int64, K::Float64, root::Box; budget=typemax(Int),
                   obbt=true, log::Union{Nothing,IO}=nothing, maxunresolved=20,
                   obbt_vars=OBBT_VARS, filter::Bool=false, excl::Vector{Box}=Box[])
    stack = [start]
    counter = Int64(0)
    nodes = 0; pruned = 0; unres = 0; maxdepth = 0; vol = 0.0
    unresolved = Tuple{Int64,Box}[]
    open_nodes = TaskNode[]
    while !isempty(stack)
        if nodes >= budget || unres >= maxunresolved
            append!(open_nodes, stack)
            break
        end
        n = pop!(stack)
        counter += 1
        id = task << 32 + counter
        res = process_node!(R, n.box, K, root; obbt, obbt_vars, filter, excl)
        nodes += 1
        maxdepth = max(maxdepth, n.depth)
        log === nothing || write_task_record(log, id, n, res)
        if res.outcome == :pruned || res.outcome == :excluded
            pruned += 1; vol += 2.0^(-n.depth)
        elseif res.outcome == :unresolved
            unres += 1; push!(unresolved, (id, res.box))
        else
            q = res.branchvar
            b1 = copy(res.box); b1.hi[q] = res.split
            b2 = copy(res.box); b2.lo[q] = res.split
            push!(stack, TaskNode(b2, n.depth + 1, id, Int8(2)))
            push!(stack, TaskNode(b1, n.depth + 1, id, Int8(1)))
        end
    end
    return (nodes=nodes, pruned=pruned, unresolved=unresolved, maxdepth=maxdepth,
            volume=vol, open=open_nodes, lps=R.lpcount, lpfail=R.lpfail)
end

# ---------------------------------------------------------------------------------------
# Depth-first search of one subtree.  `log` (optional) receives one record per node.
# ---------------------------------------------------------------------------------------
mutable struct SearchStats
    nodes::Int
    pruned::Int
    unresolved::Int
    maxdepth::Int
    volume_done::Float64
end
SearchStats() = SearchStats(0, 0, 0, 0, 0.0)

function write_record(io, id, parent, depth, res::NodeResult)
    write(io, Int64(id), Int64(parent), Int32(depth))
    code = res.outcome == :pruned ? Int8(0) : res.outcome == :branch ? Int8(1) :
           res.outcome == :excluded ? Int8(3) : Int8(2)
    write(io, code, Int32(res.branchvar), res.split)
    for v in 1:NBOX
        write(io, res.box.lo[v], res.box.hi[v])
    end
end

function search_dfs!(R::Relax, start::Box, K::Float64, root::Box; maxnodes=typemax(Int),
                     obbt_vars=OBBT_VARS, filter::Bool=false, excl::Vector{Box}=Box[],
                     obbt=true, log::Union{Nothing,IO}=nothing, stats=SearchStats(),
                     unresolved_boxes=Box[], depth0=0, idbase=0, report=0, maxunresolved=100)
    stack = [(start, depth0, Int64(-1))]
    nextid = Int64(idbase)
    open_left = Tuple{Box,Int}[]
    t0 = time()
    while !isempty(stack)
        if stats.nodes >= maxnodes
            for (b, d, _) in stack
                push!(open_left, (b, d))
            end
            break
        end
        box, depth, parent = pop!(stack)
        stats.nodes += 1
        nextid += 1
        id = nextid
        res = process_node!(R, box, K, root; obbt, obbt_vars, filter, excl)
        stats.maxdepth = max(stats.maxdepth, depth)
        log === nothing || write_record(log, id, parent, depth, res)
        if res.outcome == :pruned || res.outcome == :excluded
            stats.pruned += 1
            stats.volume_done += 2.0^(-depth)
        elseif res.outcome == :unresolved
            stats.unresolved += 1
            push!(unresolved_boxes, res.box)
            stats.unresolved >= maxunresolved && break
        else
            q = res.branchvar
            b1 = copy(res.box); b1.hi[q] = res.split
            b2 = copy(res.box); b2.lo[q] = res.split
            push!(stack, (b2, depth + 1, id))
            push!(stack, (b1, depth + 1, id))
        end
        if report > 0 && stats.nodes % report == 0
            @printf("  nodes=%d stack=%d depth<=%d done=%.3e lps=%d fail=%d %.1fs\n", stats.nodes,
                    length(stack), stats.maxdepth, stats.volume_done, R.lpcount, R.lpfail, time() - t0)
            flush(stdout)
        end
    end
    return stats, open_left, unresolved_boxes
end

end # module
