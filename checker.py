# Uses checklist.md as guide
# Any numbered checks refer to items in "Well-formed Tree"

import struct, glob, os, math, json, sys
from fractions import Fraction
from types import SimpleNamespace

K = Fraction(sys.argv[1])
LOGDIR, CERT = sys.argv[2], sys.argv[3]
KF = float(K)
assert Fraction(KF) == K
LOGS = sorted(glob.glob(os.path.join(LOGDIR, "*.bin")))

REC = struct.Struct("<qqbibid48d")
assert REC.size == 418
ALLOW_CONSTRICT = [11, 12, 13, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]

def floor_star(x):
    return Fraction(math.floor(math.nextafter(x, -math.inf) * 2**30), 2**30)

## 9. ROOT BOX exact
B0_LO = ([Fraction(0)]*4 + [Fraction(-1)]*6 + [floor_star(KF/4), floor_star(KF/2.25), floor_star(KF/2), K] + [Fraction(-1)]*4 + [Fraction(-2)]*3 + [Fraction(-9,4)]*2 + [Fraction(-4)])
B0_HI = ([Fraction(1)]*10 + [Fraction(2), Fraction(9, 4), Fraction(4), Fraction(81, 16)] + [Fraction(1)]*4 + [Fraction(2)]*3 + [Fraction(9,4)]*2 + [Fraction(4)])
assert B0_LO[10] <= K/4 and B0_LO[11] <= 4*K/9 and B0_LO[12] <= K/2, "root box invalid"

def read_records(path):
    out = []
    with open(path, "rb") as f:
        data = f.read()
    for off in range(0, len(data), REC.size):
        fields = REC.unpack_from(data, off)
        id, parent, side, depth, code, q, split = fields[:7]
        box = fields[7:] # consider 48 doubles of lo1 hi1 
        lo = [Fraction(box[2*v]) for v in range(24)]
        hi = [Fraction(box[2*v+1]) for v in range(24)]
        out.append(SimpleNamespace(id=id, parent=parent, side=side, depth=depth,
                        code=code,q=q, split=Fraction(split), 
                        lo=lo, hi=hi))
    return out

def load_tree():
    recs = {}
    for path in LOGS:
        for r in read_records(path):
            assert r.id not in recs, "duplicate id"
            recs[r.id] = r
    return recs

tree = load_tree()
print(f"N = {len(tree)}")

## 1. Root Check
root = [r for r in tree.values() if r.parent == -1]
assert len(root) == 1
assert root[0].depth == 0
assert root[0].side == 0
root = root[0]

## 2. Basic check for all other records
depthExpSum = 0
for node_id in tree:
    node = tree[node_id]
    if node.parent == -1:
        continue
    assert node.parent in tree, "parent not in tree"
    assert node.depth == tree[node.parent].depth+1, "mismatching depths"
    assert node.side == 1 or node.side == 2
    if node.code == 0:
        depthExpSum += Fraction(1, 2) ** node.depth
        assert node.q == 0, "leaf must have q = 0"
## 4c)
assert depthExpSum == 1

## Creating Children
children = {}
for node in tree.values():
    if (node.parent == -1):
        continue
    if node.parent not in children:
        children[node.parent] = [node.id]
    else:
        children[node.parent].append(node.id)

## 4a) and 4b)
vis = set()
def traverse(nid):
    assert nid not in vis, "cycle"
    vis.add(nid)
    for x in children.get(nid, []):
        traverse(x)
traverse(root.id)
assert len(tree) == len(vis), "not every node reachable from root"

## returns the prebox of 24 dimensions according to parent
def prebox(node):
    if node is root:
        return list(B0_LO), list(B0_HI)
    pa = tree[node.parent]
    lo = list(pa.lo)
    hi = list(pa.hi)
    if node.side == 1:
        hi[pa.q-1] = pa.split
    else:
        lo[pa.q-1] = pa.split
    return lo, hi

needed = []
for node_id in tree:
    ## 3. limitations on .code
    node = tree[node_id]
    assert node.code == 1 or node.code == 0, "invalid code"
    if node.code == 0:
        assert node_id not in children, "wrong for code 0"
    else: 
        assert node_id in children, "wrong for code 1"
        assert len(children[node_id]) == 2, "code 1 needs two children"
        assert sorted(tree[c].side for c in children[node_id]) == [1, 2], "only upper and lower allowed"

    ## 5. Split Validity
    if node.code == 1:
        q = node.q
        assert 1 <= q and q <= 24, "impossible q choice"
        assert node.lo[q-1] < node.split and node.split < node.hi[q-1]

    ## 6. Pre-box containment
    ## 7. ALLOW_CONSTRICT only
    pre_box_lo, pre_box_hi = prebox(node)
    for i in range(24):
        assert pre_box_lo[i] <= node.lo[i]
        assert pre_box_hi[i] >= node.hi[i]
        if (pre_box_lo[i] < node.lo[i] or pre_box_hi[i] > node.hi[i]):
            assert i+1 in ALLOW_CONSTRICT
            ## 7. append to check if all shrinking exist later
            if pre_box_lo[i] < node.lo[i]:
                needed.append([node_id, i+1, "lo"])
            if pre_box_hi[i] > node.hi[i]:
                needed.append([node_id, i+1, "hi"])

print(f"leaves = {depthExpSum == 1}, tightenings needed = {len(needed)}")

# WELL-FORMED TREE CHECKS COMPLETE! (except existence for 8&9)

def L(i, t):
    assert 5 >= i and i > t and t >= 1
    if (t == 1): return i-1
    if (t == 2): return i+2
    if (t == 3): return i+4
    if (t == 4): return i+5

def P(t):
    assert 2 <= t and t <= 5
    return t+9

def U(t, j):
    assert 1 <= t and t <= j and j <= 5 and (t, j) != (1, 1)
    if (j == t): return P(t)
    if (t == 1): return j+13
    if (t == 2): return j+16
    if (t == 3): return j+18
    if (t == 4): return j+19

TRIPLES = [(i, t, j) for t in range(1, 5) for i in range(t+1, 6) for j in range(t, 6) if (t, j) != (1, 1)]
WIDX = {trip: 25 + k for k, trip in enumerate(TRIPLES)}
def W(i, t, j): 
    return WIDX[(i, t, j)]

def TAU(i, t, j):
    if (t == i): return U(i,j)
    else: return W(i, t, j)

def var_bounds(lo24, hi24):
    lo60 = list(lo24)
    hi60 = list(hi24)
    for trip in TRIPLES:
        i, t, j = trip
        a = lo24[L(i, t)-1]
        b = hi24[L(i, t)-1]
        c = lo24[U(t, j)-1]
        d = hi24[U(t, j)-1]
        mn = min({a*c, a*d, b*c, b*d})
        mx = max({a*c, a*d, b*c, b*d})
        lo60.append(mn)
        hi60.append(mx)
    assert len(lo60) == 60 and len(hi60) == 60
    for m in range(60): assert lo60[m] <= hi60[m]
    return lo60, hi60

CONST_R = [0.5, 0.75, 1.0, 1.25, 1.375, 1.5, 1.625, 1.75, 2.0]

def build_rows(lo24, hi24):
    rows = []
    def add(name, coefs, rhs):
        rows.append((name, {k: Fraction(v) for k, v in coefs.items()}, Fraction(rhs)))

    # CP
    for i in range(1, 6):
        for j in range(2, 6):
            coefsneg = {}
            coefspos = {}
            for t in range(1, min(i, j)+1):
                coefsneg[TAU(i, t, j)] = -1
                coefspos[TAU(i, t, j)] = 1
            rhs = 1
            add(f"CP({1},{i},{j},-)", coefsneg, rhs)
            add(f"CP({1},{i},{j},+)", coefspos, rhs)

    for k in range(2, 5):
        for i in range(k, 6):
            for j in range(k+1, 6):
                coefsneg = {}
                coefspos = {}
                for t in range(k, min(i, j)+1):
                    coefsneg[TAU(i, t, j)] = -1
                    coefspos[TAU(i, t, j)] = 1
                coefsneg[P(k)] = -1
                coefspos[P(k)] = -1
                add(f"CP({k},{i},{j},-)", coefsneg, 0)
                add(f"CP({k},{i},{j},+)", coefspos, 0)

    # COL
    for t in range(2, 5):
        for i in range(t+1, 6):
            add(f"COL({i},{t},-)", {P(t): -1, W(i, t, t): -1}, 0)
            add(f"COL({i},{t},+)", {P(t): -1, W(i, t, t): 1}, 0)

    # PIVOT: DBL(k)
    add("DBL(1)", {P(2): 1}, 2)
    for k in range(2, 5):
        add(f"DBL({k})", {P(k+1): 1, P(k): -2}, 0)

    # PIVOT: COHEN(k, r)
    for k in range(1, 4):
        for r in CONST_R:
            # p{k+2} - r^2 p_k - (3-2r) p_{k+1} <= 0
            if (k == 1): add(f"COHEN(1,{r})", {P(3): 1, P(2): 2*r-3}, r*r)
            else: add(f"COHEN({k},{r})", {P(k): -r**2, P(k+1): 2*r-3, P(k+2): 1}, 0)

    # PIVOT: G4
    add(f"G4", {P(5):1, P(2):-4}, 0)

    # McCormick rows: MC(i,t,j,lolo/lohi/hilo/hihi)
    for trip in TRIPLES:
        i, t, j = trip
        a = lo24[L(i, t)-1]
        b = hi24[L(i, t)-1]
        c = lo24[U(t, j)-1]
        d = hi24[U(t, j)-1]
        l = L(i, t)
        u = U(t, j)

        add(f"MC({i},{t},{j},lolo)", {l: c, u: a, W(i, t, j): -1}, a*c)
        add(f"MC({i},{t},{j},lohi)", {l: -d, u: -a, W(i, t, j): 1}, -a*d)
        add(f"MC({i},{t},{j},hilo)", {l: -c, u: -b, W(i, t, j): 1}, -b*c)
        add(f"MC({i},{t},{j},hihi)", {l: d, u: b, W(i, t, j): -1}, b*d)

    # Perspective rows: PERSP(i,t,j,pm)
    for trip in TRIPLES:
        i, t, j = trip
        if j <= t: continue
        a = lo24[L(i, t)-1]
        b = hi24[L(i, t)-1]
        c = lo24[U(t, j)-1]
        d = hi24[U(t, j)-1]
        u = U(t, j)
        w = W(i, t, j)
        m = (a+b)/2
        r = (b-a)/2
        if t == 1:
            add(f"PERSP({i},{t},{j},-)", {w:-1, u:m}, r)
            add(f"PERSP({i},{t},{j},+)", {w:1, u:-m}, r)
        else:
            add(f"PERSP({i},{t},{j},-)", {w:-1, u:m, P(t):-r}, 0)
            add(f"PERSP({i},{t},{j},+)", {w:1, u:-m, P(t):-r}, 0)

    assert len(rows) == 328, "wrong number of rows"
    assert len({r[0] for r in rows}) == 328, "duplicate row name"
    return rows

rows0 = build_rows(B0_LO, B0_HI)
lo60, hi60 = var_bounds(B0_LO, B0_HI)
print(f"rows = {len(rows0)}, w212 = [{lo60[W(2,1,2)-1]}, {hi60[W(2,1,2)-1]}], w322 = [{lo60[W(3,2,2)-1]}, {hi60[W(3,2,2)-1]}]")

# rows, list, list dict, list;
def Q(rows, lo60, hi60, c, lam):
    d = [Fraction(0)]*60
    for var, coef in c.items():
        d[var-1] = coef
    htl = Fraction(0)
    for name, coefs, rhs in rows:
        lam_i = lam.get(name, 0)
        if lam_i == 0: continue
        htl += lam_i * rhs
        for var, coef in coefs.items():
            d[var-1] -= lam_i * coef

    suma = 0
    for i in range(60):
        suma += max(d[i]*lo60[i], d[i]*hi60[i])
    return htl + suma

VARIDX = {}
for t in range(1, 5):
    for i in range(t+1, 6):
        VARIDX[f"L{i}{t}"] = L(i, t)
for t in range(2, 6):
    VARIDX[f"p{t}"] = P(t)
for t in range(1, 5):
    for j in range(t+1, 6):
        VARIDX[f"U{t}{j}"] = U(t, j)
assert len(VARIDX) == 24 and sorted(VARIDX.values()) == list(range(1, 25))

needed_set = {(nid, v, s) for nid, v, s in needed}
certified_leaves = set()
cache_id, cache = None, None
cnt = 0

def node_data(nid):
    global cache_id, cache
    if nid != cache_id:
        node = tree[nid]
        pre_lo, pre_hi = prebox(node)
        cache_id = nid
        cache = (node, build_rows(pre_lo, pre_hi), *var_bounds(pre_lo, pre_hi),
                 build_rows(node.lo, node.hi), *var_bounds(node.lo, node.hi))
    return cache

# All exact calculations for shrinking and pruning are made here, 
# "Certificates" section in checklist.md
with open(CERT) as f:
    header = json.loads(next(f))
    for line in f:
        if cnt % 100000 == 0:
            print(f"{cnt} Certificates Checked")
        cnt += 1
        cert = json.loads(line)
        nid = cert["node"]
        assert nid in tree, f"certificate for unknown node {nid}"
        node, pre_rows, pre_lo60, pre_hi60, box_rows, box_lo60, box_hi60 = node_data(nid)
        lam = {name: Fraction(x) for name, x in cert["lam"]}
        assert all(x > 0 for x in lam.values())
        if cert["kind"] == "leaf":
            assert node.code == 0
            Qc0 = Q(box_rows, box_lo60, box_hi60, {}, lam)
            assert Qc0 == Fraction(cert["Q"]), f"Q mismatch at {nid}"
            assert Qc0 < 0, f"leaf not certified at {nid}"
            certified_leaves.add(nid)
        else:
            v = VARIDX[cert["var"]]
            side = cert["side"]
            coef = Fraction(cert["objective"]["coef"])
            assert coef == (1 if side == "hi" else -1), "wrong side"
            Qcl = Q(pre_rows, pre_lo60, pre_hi60, {v: coef}, lam)
            assert Qcl == Fraction(cert["Q"]), f"Q mismatch at {nid} {cert['var']} {side}"
            target = node.hi[v-1] if side == "hi" else -node.lo[v-1]
            assert target == Fraction(cert["target"])
            assert Qcl <= target, "invalid shrinking"
            needed_set.discard((nid, v, side))

leaf_ids = {nid for nid, n in tree.items() if n.code == 0}
## remaining parts of 8. and 9.
assert leaf_ids == certified_leaves, "leaf without certificate"
assert not needed_set, "tightening without certificate"
print(f"verified {len(needed)} tightenings and {len(certified_leaves)} leaves")
