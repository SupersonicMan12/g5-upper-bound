# Data formats

This file gives the layouts of the files that `checker.py` reads: the node logs in
`data/logs_K*/` and the certificate files `data/certificates_K*.jsonl`. It describes
formats only; the mathematics (the 24 coordinates, the relaxation and its rows) is in the
paper.

## Coordinates

Boxes have 24 coordinates, numbered 1–24 in this order:

| index | coordinates |
|---|---|
| 1–10 | $\ell_{21},\ell_{31},\ell_{41},\ell_{51},\ \ell_{32},\ell_{42},\ell_{52},\ \ell_{43},\ell_{53},\ \ell_{54}$ |
| 11–14 | $p_2,p_3,p_4,p_5$ |
| 15–24 | $u_{12},u_{13},u_{14},u_{15},\ u_{23},u_{24},u_{25},\ u_{34},u_{35},\ u_{45}$ |

In certificates they are named `L21`, ..., `L54`, `p2`, ..., `p5`, `U12`, ..., `U45`. The
relaxation adds 36 product variables $w_{itj}$ (indices 25–60), for $1\le t\le4$,
$t<i\le5$, $t\le j\le5$, $(t,j)\ne(1,1)$, ordered by $t$, then $i$, then $j$.

## Node logs: `data/logs_K*/task_<n>.bin`

Each file is a sequence of fixed-size records, 418 bytes each, little-endian. In Python:
`struct.Struct("<qqbibid48d")`.

| offset | bytes | type | field |
|---|---|---|---|
| 0 | 8 | int64 | `id` = `task << 32` + counter (unique over the whole run) |
| 8 | 8 | int64 | `parent` id (−1 for the root) |
| 16 | 1 | int8 | `side`: 0 root; 1 lower child (`hi[q] = split`); 2 upper child (`lo[q] = split`) |
| 17 | 4 | int32 | `depth` (root = 0) |
| 21 | 1 | int8 | `code`: 0 pruned leaf; 1 branch node |
| 22 | 4 | int32 | `q`: branch coordinate (1–24) if `code = 1`; else 0 |
| 26 | 8 | float64 | `split`: split point (meaningful only if `code = 1`) |
| 34 | 384 | 24 × (float64, float64) | the node's box **after tightening**: `(lo, hi)` per coordinate, in the order above |

The search code also defines codes 2 (unresolved leaf) and 3 (excluded leaf); they do not
occur in these logs, and `checker.py` rejects them.

**Semantics.**

- **Boxes before tightening.** The *pre-box* of the root is the root box $B_0(K)$, which is
  not stored. The pre-box of any other node is its parent's stored box with `hi[q] = split`
  (side 1) or `lo[q] = split` (side 2), where `q` and `split` are the parent's.
- **Tightening.** The stored box must lie inside the pre-box. Only $p_2,p_3,p_4$ and the
  ten $u$'s (coordinates 11–13 and 15–24) may be tightened.
- **Leaves.** `code = 0`: no point of the relaxation of the stored box has $p_5\ge K$. This
  is what a leaf certificate shows.
- **Order.** Records of one task are in depth-first order. The first record of a task may
  have its parent in another file, because large tasks are split.

Sizes: `data/logs_K4p5` holds 7,543 records in 73 files; `data/logs_K4p15` holds 225,465
records in 1,019 files.

## The root box $B_0(K)$, as computed by the code

| coordinates | interval |
|---|---|
| $\ell_{i1}$ | $[0,1]$ |
| $\ell_{it}$, $t\ge2$ | $[-1,1]$ |
| $p_2,\ p_3,\ p_4,\ p_5$ | $[\lfloor K/4\rfloor_*,\,2]$, $[\lfloor K/2.25\rfloor_*,\,2.25]$, $[\lfloor K/2\rfloor_*,\,4]$, $[K,\,5.0625]$ |
| $u_{1j}$ | $[-1,1]$ |
| $u_{tj}$, $t\ge2$ | $[-\bar p_t,\bar p_t]$ with $\bar p_2=2$, $\bar p_3=2.25$, $\bar p_4=4$ |

Here $\lfloor x\rfloor_*$ = `floor(prevfloat(fl(x)) * 2^30) / 2^30`, evaluated in double
precision. `fl(x)` is the computed double quotient, and `prevfloat` is the next smaller
double. A checker must reproduce these doubles exactly, so that they match the logs, and
must prove that each lower bound is at most the exact value $K/4$, $4K/9$ or $K/2$.

**$K$.** The search code rounds the requested value down to a multiple of $2^{-20}$:
$K = 9/2$ for 4.5 and $K = 2175795/2^{19}$ for 4.15. The search reports give both
`K_dyadic` and `K_rational`.

## Certificate files: `data/certificates_K4p5.jsonl`, `data/certificates_K4p15.jsonl`

Produced by `generators/export_v3/lubb_export.jl`, with 0 failures:

| file | size | certificates | tightenings | leaves |
|---|---|---|---|---|
| `certificates_K4p5.jsonl` | 154 MB | 102,899 | 99,127 | 3,772 |
| `certificates_K4p15.jsonl` | 6.18 GB | 3,478,969 | 3,366,236 | 112,733 |

JSON lines: the first line is a header, every other line is one certificate. Read the file
with `json.loads(line)` per line (do not load the whole file into one JSON array). Every
multiplier was proposed by HiGHS and then checked exactly in dyadic arithmetic before being
written; the files contain **only** certificates that passed.

**Header** (line 1): `format`, `K` (as a double), `K_rational`, `logdir` (the log path at
export time, e.g. `4.5/logs_K4p5`; same data as `data/logs_K4p5`), `nodes`, `tightenings`,
`leaves`, `failures` (0), `lubb_sha256` (the SHA-256 of `generators/export_v3/lubb.jl`),
`date_utc`.

**Certificate** (every other line):

| key | meaning |
|---|---|
| `node` | the node `id` of the log record |
| `kind` | `"tighten"` or `"leaf"` |
| `var`, `side`, `target` | (`tighten` only) the coordinate (`L21`, `p3`, `U24`, ...), `"hi"` or `"lo"`, and the certified value: `target` $=$ the stored `hi[v]` for `hi`, $= -$`lo[v]` for `lo` |
| `objective` | `null` ($c=0$) or `{"var": name, "coef": ±1.0}` ($c=\pm e_v$) |
| `check` | what the exact value $Q$ of the LP-duality bound (see the paper) must satisfy: `"Q<0"` or `"Q<=target"` |
| `Q` | the exact value as `"numerator/denominator"`, computed by the exporter; an independent computation must agree **exactly** |
| `lam` | list of `[rowname, multiplier]` for the rows with multiplier $>0$; all other rows have multiplier $0$ |

**Which box.** A `tighten` certificate is evaluated on the node's **pre-box** (the parent's
stored box split at the logged point; $B_0(K)$ for the root). A `leaf` certificate is
evaluated on the node's **stored** (tightened) box.

**Numbers.** Multipliers and `target` are doubles printed as the shortest decimal that
reads back to the same double; `Fraction(x)` of the parsed float is the exact value. Every
leaf has `objective = null` and `check = "Q<0"`, and every tightening has
`check = "Q<=target"`. The multipliers come from version 3 of `lubb.jl` with the version-1
options (theorems on, no symmetry rows), not from the search code; the exact check is what
makes them valid.

**Row names.** Rows are identified by name, not position, so a checker may order them as it
likes. Each name denotes one inequality, written here as `expression <= rhs` with
coefficients exactly as shown (a multiplier only matches the row at this scaling).
$\tau_{itj}$ means $u_{ij}$ when $t=i$ and $w_{itj}$ otherwise; the box gives
$\ell_{it}\in[a,b]$, $u_{tj}\in[c,d]$, $m=(a+b)/2$, $r=(b-a)/2$.

| name | inequality | count |
|---|---|---|
| `CP(k,i,j,±)` for $1\le k\le4$, $k\le i\le5$, $k<j\le5$ | $\pm\sum_{t=k}^{\min(i,j)}\tau_{itj} - p_k\le0$; for $k=1$: $\pm\sum\le1$ | 80 |
| `COL(i,t,±)` for $2\le t<i\le5$ | $\pm w_{itt}-p_t\le0$ | 12 |
| `DBL(k)` for $1\le k\le4$ | $p_{k+1}-2p_k\le0$; for $k=1$: $p_2\le2$ | 4 |
| `COHEN(k,r)` for $1\le k\le3$, $r\in\{0.5,0.75,1.0,1.25,1.375,1.5,1.625,1.75,2.0\}$ | $p_{k+2}-r^2p_k-(3-2r)p_{k+1}\le0$; for $k=1$: $p_3-(3-2r)p_2\le r^2$ | 27 |
| `G4` | $p_5-4p_2\le0$ | 1 |
| `MC(i,t,j,lolo)` | $-w+c\,\ell+a\,u\le ac$ | 36 |
| `MC(i,t,j,hihi)` | $-w+d\,\ell+b\,u\le bd$ | 36 |
| `MC(i,t,j,lohi)` | $w-d\,\ell-a\,u\le-ad$ | 36 |
| `MC(i,t,j,hilo)` | $w-c\,\ell-b\,u\le-bc$ | 36 |
| `PERSP(i,t,j,±)` for $j>t$ | $\pm(w-m\,u)-r\,p_t\le0$; for $t=1$: $\pm(w-m\,u)\le r$ | 60 |

Here $w=w_{itj}$, $\ell=\ell_{it}$, $u=u_{tj}$ (which is $p_t$ when $j=t$). `r` in
`COHEN` is printed as a decimal (`1.375`), and `±` as `+` or `-`. Total: 328 rows.

**Checking a certificate.** Build the 328 rows and the 60 variable bounds from the box;
put the listed multipliers in your row order (zero elsewhere); evaluate the exact bound
with the stated objective; compare with `Q`; then apply `check`.
