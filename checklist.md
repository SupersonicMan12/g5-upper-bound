# Well-formed Tree

Node $\nu$ has pre-box $B_\nu = [\underline b, \overline b]$ and stored box $B'_\nu = [\underline b', \overline b']$, with coordinates $v = 1, \dots, 24$.

1. Exactly one record with parent = -1, depth = 0, side = 0

2. All other records have proper `parent_id`, side either 1 or 2, depth = parent.depth + 1

3. Code-1 nodes have exactly two children with sides 1 and 2, code-0 should not have any children

4. Tree properties: a) no cycles; b) every node reachable from root; c) Exactly `sum(Fraction(1, 2**r.depth) for leaves) == 1`, although redundant given the previous

5. Split validity: $q$ within 1 and 24 and $\underline b'_q < \text{split} < \overline b'_q$; if equal to one side technically no problem, but will be flagged

6. Pre-box containment, ensure that $B'_\nu \subseteq B_\nu$

7. Whenever $\underline b'_v > \underline b_v$ or $\overline b'_v < \overline b_v$ there is $v \in \{11, 12, 13, 15, \cdots, 24\}$ and record `(id, v, side)` so we can check if the corresponding shrinking certificate exists

8. At any leaf, check that the corresponding pruning certificate exists

9. Root box is exactly $B_0(K)$, with lower bounds of $p_2, p_3, p_4$ computed by $\lfloor\cdot\rfloor_*$, and these are at most $K/4$, $4K/9$, $K/2$

# Certificates

Node $\nu$ with pre-box $B$ and stored box $B'$. Let $R(B)$ be the 60-dimension relaxation and the $Gx \le h$. Define for $\lambda \ge 0$ the objective $c$:

$$
Q(c, \lambda) = \lambda^T h + \sum_{v=1}^{60} \max(d_v\underline x_v, d_v \overline x_v)
$$
for $d = c-G^T\lambda$. We have that $c^T x \le Q(c, \lambda)$ for every $x\in R(B)$.

## Shrinking is Valid

On $R(B)$:

$$
Q(e_v, \lambda) \le \overline b'_v\\
Q(-e_v, \lambda) \le -\underline b'_v
$$
This ensures $\mathcal F_K \cap B \subseteq B'$.

## Pruning is Valid

$Q(0,\lambda) < 0$ on $R(B')$ (the stored box).
