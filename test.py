
import sympy

def sort_roots(roots):
    roots = [(root.as_real_imag()) for root in roots]
    real_roots = [root for root in roots if root[1] == 0]
    real_roots = sorted(real_roots)
    com_roots = [root for root in roots if root[1] != 0]
    com_roots  = sorted(com_roots) ## does sympy have sorted :|
    return real_roots + com_roots


def make_algebraic_fsymbol(string):
    coefficients = [int(x) for x in string[:-3].split("_")]
    index = int(string[-1])-1 

    x = sympy.Symbol('x')
    expression = sum(c * x**i for i, c in enumerate(coefficients))
    

print(make_algebraic_fsymbol("1_1_1__2"))
