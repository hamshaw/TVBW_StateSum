"""
Turaev-Viro-Barrett-Westbury state sum for a closed 3-manifold triangulation
(Regina) and a spherical fusion category (JSON file, looked up by UUID).

Conventions (multiplicity-free or not), unitary gauge assumed:
  F key  "a_b_c_d_e_mu_nu_f_rho_sigma"  (element "index"), meaning
      (a x b -> e; mu) , (e x c -> d; nu)   [left  bracketing]
      (b x c -> f; rho), (a x f -> d; sigma) [right bracketing]
  Object 1 is the unit.
  In multiplicity-free case, mu, nu, rho, sigma = 0

State sum (vertices ordered by a branching structure i<j<k<l in each tet):
  Z = D^{-2V} * sum_{edge colours, triangle multiplicity labels}
        prod_edges d_e * prod_tets  F^{a_ij a_jk a_kl}_{a_il; a_ik, a_jl} / sqrt(d_ik d_jl)
  (complex conjugate for tets whose branching orientation disagrees with the
   manifold orientation),   D^2 = sum_a d_a^2.
With this normalisation Z(S^3) = 1/dim(C), Z(S^2 x S^1) = 1.

All category data (F-symbols, quantum dimensions) is kept as EXACT algebraic
numbers via sympy, parsed from the (non-"numeric_") string fields of the
AnyonWiki JSON format, rather than as floating-point complex numbers.

AnyonWiki algebraic-number string format (see
https://github.com/anyonwiki/AnyonWikiDatabase):
    "c0_c1_..._cn__k"
  c0..cn are the integer coefficients (constant term first) of a minimal
  polynomial c0 + c1*x + ... + cn*x^n, and k (after the "__") is the 1-based
  index into that polynomial's roots taken in a fixed canonical order (real
  roots ascending, then complex roots ordered by ascending real part, then
  ascending imaginary part -- the same order sympy's Poly.all_roots() uses).
  Example: "-1_4__1" -> -1 + 4x = 0 -> x = 1/4.
           "1_0_1__2" -> x^2 + 1 = 0, roots [-I, I] -> pick #2 -> I.
"""
import json
import itertools
from collections import defaultdict

import sympy as sp
import regina

_x = sp.symbols("x")
_root_cache = {}


def parse_algebraic(s): # Parse one AnyonWiki algebraic-number string into an exact sympy number.
    if s in _root_cache:
        return _root_cache[s] #optimize finding roots
    poly_part, idx_part = s.rsplit("__", 1)
    idx = int(idx_part)                       # 1-based
    coeffs = [int(c) for c in poly_part.split("_")]       # c0 (const) .. cn
    poly = sp.Poly(list(reversed(coeffs)), _x)            # sympy wants leading-first
    roots = poly.all_roots()                  # canonical, deterministic order
    val = sp.simplify(roots[idx - 1])
    _root_cache[s] = val
    return val


# ----------------------------------------------------------------------------
# Category
#   Store necessary data from the SFC
# ----------------------------------------------------------------------------
class SphericalFusionCategory:
    def __init__(self, cat_id, json_path, fusion_rules=None):
        """
        fusion_rules: optional dict {(a,b,c): N_ab^c} with 0-based objects.
        If None, N is inferred from the F-symbol keys (fine as long as the
        JSON lists F-symbols involving every admissible channel).
        Use None for now
        """
        with open(json_path) as fh:
            data = json.load(fh)
        entry = data["data"][cat_id] if cat_id in data["data"] else data
        self.id = cat_id

        self.dims = [parse_algebraic(s) for s in entry["quantum_dimensions"]]
        self.rank = len(self.dims)
        self.dim = sp.nsimplify(sum(d * d for d in self.dims))   # D^2 = dim(C)
        self.is_spherical = entry["is_spherical"]
        self.F = {}
        for key, val in entry["f_symbols"].items():
            a, b, c, d, e, mu, nu, f, rho, sg = [int(x) - 1 for x in key.split("_")]
            self.F[(a, b, c, d, e, f, mu, nu, rho, sg)] = parse_algebraic(val)

        if fusion_rules is not None:
            self.N = defaultdict(int, fusion_rules)
        else:
            self.N = defaultdict(int)
            for (a, b, c, d, e, f, mu, nu, rho, sg) in self.F:
                self.N[(a, b, e)] = max(self.N[(a, b, e)], mu + 1)
                self.N[(e, c, d)] = max(self.N[(e, c, d)], nu + 1)
                self.N[(b, c, f)] = max(self.N[(b, c, f)], rho + 1)
                self.N[(a, f, d)] = max(self.N[(a, f, d)], sg + 1)

        self._sqrt_cache = {}

    def _inv_sqrt(self, e, f): # to calculate the leading normalization factor
        key = (e, f) if e <= f else (f, e)
        val = self._sqrt_cache.get(key)
        if val is None:
            val = 1 / sp.sqrt(self.dims[e] * self.dims[f])
            self._sqrt_cache[key] = val
        return val

    def tet_weight(self, a, b, c, d, e, f, mu, nu, rho, sg): # gets tetrahedra weight given F key
        v = self.F.get((a, b, c, d, e, f, mu, nu, rho, sg), sp.Integer(0))
        if v == 0:
            return sp.Integer(0)
        return v * self._inv_sqrt(e, f)


# ----------------------------------------------------------------------------
# Branching structure on a Regina triangulation
# ----------------------------------------------------------------------------
def _tet_local_edges(tri):
    """For each tet: list of (regina_edge_index, u, w) with (u,w) = edge vertices
    in the edge's own vertex order (vertex 0 -> vertex 1), as tet-local labels."""
    out = []
    for t in tri.simplices():
        lst = []
        for k in range(6):
            m = t.edgeMapping(k)
            lst.append((t.edge(k).index(), m[0], m[1]))
        out.append(lst)
    return out


def find_branching(tri):
    """Orient each edge so that every tet gets a total order on its vertices.
    Returns list s[e] in {0,1} (0: edge points from its vertex 0 to vertex 1)
    or None if none exists."""
    tel = _tet_local_edges(tri)
    ne = tri.countEdges()
    # order edges by first appearance so tets close early
    order, seen = [], set()
    for lst in tel:
        for (e, _, _) in lst:
            if e not in seen:
                seen.add(e)
                order.append(e)
    pos = {e: i for i, e in enumerate(order)}
    closes = defaultdict(list)
    for ti, lst in enumerate(tel):
        closes[max(pos[e] for (e, _, _) in lst)].append(ti)

    s = [None] * ne

    def acyclic(ti):
        outdeg = [0] * 4
        for (e, u, w) in tel[ti]:
            tail = u if s[e] == 0 else w
            outdeg[tail] += 1
        return sorted(outdeg) == [0, 1, 2, 3]

    def rec(i):
        if i == ne:
            return True
        e = order[i]
        for val in (0, 1):
            s[e] = val
            if all(acyclic(t) for t in closes[i]) and rec(i + 1):
                return True
        s[e] = None
        return False

    return s if rec(0) else None


def branched_data(tri):
    """Returns (tri_used, per-tet info). Falls back to a barycentric
    subdivision if no branching exists on tri."""
    s = find_branching(tri)
    if s is None:
        tri = regina.Triangulation3(tri)
        tri.subdivide()          # barycentric subdivision (always branchable)
        tri.orient()
        s = find_branching(tri)
        assert s is not None
    tel = _tet_local_edges(tri)
    info = []
    for ti, t in enumerate(tri.simplices()):
        outdeg = [0] * 4
        for (e, u, w) in tel[ti]:
            outdeg[u if s[e] == 0 else w] += 1
        p = sorted(range(4), key=lambda v: -outdeg[v])   # p[0]<p[1]<p[2]<p[3]
        # sign of permutation p
        inv = sum(1 for i in range(4) for j in range(i + 1, 4) if p[i] > p[j])
        sign = t.orientation() * (-1 if inv % 2 else 1)

        E = lambda i, j: t.edge(regina.Edge3.edgeNumber[p[i]][p[j]]).index()
        T = lambda opp: t.triangle(p[opp]).index()
        info.append(dict(
            a=E(0, 1), b=E(1, 2), c=E(2, 3), d=E(0, 3), e=E(0, 2), f=E(1, 3),
            mu=T(3), nu=T(1), rho=T(0), sigma=T(2), sign=sign))
    return tri, info


# ----------------------------------------------------------------------------
# State sum
# ----------------------------------------------------------------------------
def turaev_viro(tri, cat, verbose=False):

    # quality checks for input
    if not tri.isValid() or not tri.isClosed():
        raise ValueError("triangulation must be valid and closed")
    if not tri.isOrientable():
        raise ValueError("triangulation must be orientable")
    if not cat.is_spherical:
        raise ValueError("category must be spherical")

    # set up input
    tri = regina.Triangulation3(tri)
    tri.orient()

    tri, info = branched_data(tri)
    V, ne = tri.countVertices(), tri.countEdges()
    R = cat.rank

    # edge ordering by first appearance; tets closed when all 6 edges coloured
    order, seen = [], set()
    for t in info: # t is a tet
        for k in "abcdef": # k is one edge in that tet
            if t[k] not in seen:
                seen.add(t[k]); order.append(t[k])
    pos = {e: i for i, e in enumerate(order)}
    closes = defaultdict(list)
    for ti, t in enumerate(info):
        closes[max(pos[t[k]] for k in "abcdef")].append(ti)

    col = [None] * ne
    total = sp.Integer(0)

    def admissible(t):
        A, B, C, D_, E_, F_ = (col[t[k]] for k in "abcdef")
        N = cat.N
        return (N[(A, B, E_)] and N[(E_, C, D_)] and
                N[(B, C, F_)] and N[(A, F_, D_)]) # confirms all vertices are valid

    def leaf(): # add weights of a addmissible labeling
        nonlocal total
        # multiplicity spaces of triangles
        tri_range = {}
        for t in info:
            A, B, C, D_, E_, F_ = (col[t[k]] for k in "abcdef")
            tri_range[t["mu"]] = cat.N[(A, B, E_)]
            tri_range[t["nu"]] = cat.N[(E_, C, D_)]
            tri_range[t["rho"]] = cat.N[(B, C, F_)]
            tri_range[t["sigma"]] = cat.N[(A, F_, D_)] # always going to be 1 in multiplicity free. This works but maybe can be optimized.
        tids = list(tri_range)
        edge_w = sp.Integer(1)
        for e in range(ne):
            edge_w *= cat.dims[col[e]]
        acc = sp.Integer(0)
        for labels in itertools.product(*(range(tri_range[x]) for x in tids)):
            lab = dict(zip(tids, labels))
            w = sp.Integer(1)
            for t in info:
                A, B, C, D_, E_, F_ = (col[t[k]] for k in "abcdef")
                x = cat.tet_weight(A, B, C, D_, E_, F_, lab[t["mu"]],
                                   lab[t["nu"]], lab[t["rho"]], lab[t["sigma"]]) # again, can optimize if garunteed to be mult-free
                if x == 0:
                    w = sp.Integer(0)
                    break
                w *= sp.conjugate(x) if t["sign"] < 0 else x
            acc += w
        total += edge_w * acc

    def rec(i): # i is how many edges we've decided thus far
        if i == ne:# end if all are decided
            leaf() # add this colorings contribution to the sum
            return
        e = order[i] # current edge we care about
        for c in range(R): # try all objects in R = cat.rank
            col[e] = c
            if all(admissible(info[ti]) for ti in closes[i]):
                rec(i + 1)
        col[e] = None # failed to find an addmissible object for this edge

    rec(0) # start with 0 decided edges
    Z = sp.simplify(total / (cat.dim ** V))
    return Z


if __name__ == "__main__":
    import sys
    #cat = SphericalFusionCategory(sys.argv[1], sys.argv[2])
    ## tri = regina.Triangulation3.fromSig(sys.argv[3]) if len(sys.argv) > 3 \
        
    tri = regina.Example3.lens(5, 1)
    #tri = regina.Example3.poincare()
    cats = ['37fa4c58-97c6-11f1-a6a6-639d9c4bace6',
            '37fa4e2e-97c6-11f1-aa03-3b817be52f72', 
            '37fa4ea8-97c6-11f1-8b25-ad9b74351a3e', 
            '37fa4f78-97c6-11f1-90cf-c3d624fdb861', 
            '37fa5126-97c6-11f1-8476-3f003633fef2', 
            '37fa5642-97c6-11f1-9458-d142162199ed', 
            '37fa5c3e-97c6-11f1-8693-91f67133af96']
    for cate in cats:
        cat = SphericalFusionCategory(cate, sys.argv[2])
        Z = turaev_viro(tri, cat)
        print(f'{sys.argv[1]} \t Rank:{cat.rank} \t Sum:{Z}')
    #print("Z (numeric)=", sp.N(Z))
