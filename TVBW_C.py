"""
Turaev-Viro-Barrett-Westbury state sum for a closed 3-manifold triangulation
(Regina) and a spherical fusion category (JSON file, looked up by UUID).

Conventions (multiplicity-free or not), unitary gauge assumed:
  F key  "a_b_c_d_e_mu_nu_f_rho_sigma"  (1-based), meaning
      (a x b -> e; mu) , (e x c -> d; nu)   [left  bracketing]
      (b x c -> f; rho), (a x f -> d; sigma) [right bracketing]
  Object 1 is the unit.

State sum (vertices ordered by a branching structure i<j<k<l in each tet):
  Z = D^{-2V} * sum_{edge colours, triangle multiplicity labels}
        prod_edges d_e * prod_tets  F^{a_ij a_jk a_kl}_{a_il; a_ik, a_jl} / sqrt(d_ik d_jl)
  (complex conjugate for tets whose branching orientation disagrees with the
   manifold orientation),   D^2 = sum_a d_a^2.
With this normalisation Z(S^3) = 1/dim(C), Z(S^2 x S^1) = 1.
"""
import json
import cmath
import itertools
from collections import defaultdict

import sympy
import regina

def make_algebraic_number(string):
    coefficients = [int(x) for x in string[:-3].split("_")] # get the polynomial coefficients

    index = int(string[-1]) - 1 # get the index of the root we want from the polynomial
    if string[-2] != "_":
        index = int(string[-2:]) - 1 # accounts for it being 2 digit if necessary

    x = sympy.Symbol('x')
    expression = sum(c * x**i for i, c in enumerate(coefficients))

    return sympy.RootOf(expression, index)

# ----------------------------------------------------------------------------
# Category
# ----------------------------------------------------------------------------
class SphericalFusionCategory:
    def __init__(self, cat_id, json_path, fusion_rules=None):
        """
        fusion_rules: optional dict {(a,b,c): N_ab^c} with 0-based objects.
        If None, N is inferred from the F-symbol keys (fine as long as the
        JSON lists F-symbols involving every admissible channel).
        """
        with open(json_path) as fh:
            data = json.load(fh)
        entry = data[cat_id]# if cat_id in data else data
        self.id = cat_id
        self.dims = [make_algebraic_number(x) for x in entry["quantum_dimensions"]]
        self.rank = len(self.dims)
        self.dim = sum(d * d for d in self.dims)          # D^2 = dim(C)
        self.is_spherical = entry["is_spherical"]
        self.F = {}
        for key, val in entry["f_symbols"].items():
            a, b, c, d, e, mu, nu, f, rho, sg = [int(x)-1 for x in key.split("_")] #used to be int(x)-1??
            self.F[(a, b, c, d, e, f, mu, nu, rho, sg)] = make_algebraic_number(val)
            print(type(self.F[(a, b, c, d, e, f, mu, nu, rho, sg)]))

        if fusion_rules is not None:
            self.N = defaultdict(int, fusion_rules)
        else:
            self.N = defaultdict(int)
            for (a, b, c, d, e, f, mu, nu, rho, sg) in self.F:
                self.N[(a, b, e)] = max(self.N[(a, b, e)], mu + 1)
                self.N[(e, c, d)] = max(self.N[(e, c, d)], nu + 1)
                self.N[(b, c, f)] = max(self.N[(b, c, f)], rho + 1)
                self.N[(a, f, d)] = max(self.N[(a, f, d)], sg + 1)

        self._sqrt = {}

    def tet_weight(self, a, b, c, d, e, f, mu, nu, rho, sg):
        v = self.F.get((a, b, c, d, e, f, mu, nu, rho, sg), 0)
        if v == 0:
            return sympy.S.Zero
        s = cmath.sqrt(self.dims[e] * self.dims[f])
        return sympy.S(v / s)
# ----------------------------------------------------------------------------
# State sum
# ----------------------------------------------------------------------------
def turaev_viro(tri, cat):

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

    return 0



if __name__ == "__main__":
    import sys
    cat = SphericalFusionCategory(sys.argv[1], sys.argv[2])
    tri = regina.Triangulation3.fromSig(sys.argv[3]) if len(sys.argv) > 3 \
        else regina.Example3.lens(5, 1)
