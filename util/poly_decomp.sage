def randvect(G, size):
    return vector(G, [G.random_element() for _ in range(size)])

def poly_decomp(p, x, t):
    """
    p: The polynomial of form `f(x) - f(t)`
    x: A vector of variables in `R`
    t: A vector with same dimension as `x`

    The polynomial will be factored by (x_i - t_i). Returns a list of coefficients of each term.
    """
    R = p.parent() # polynomial ring containing p
    result = []
    for i in range(len(x)):
        d = eval(x[i]) - t[i]
        J = R.ideal(d)

        r = p.reduce(J)
        q = (p - r) / d

        result.append(q)
        p = r

    return result

"""
from Crypto.Util.number import getPrime

N = 5
q = getPrime(2048)
G = Zmod(q).field() # since q is prime, Zmod(q) is a field

variables = ["x%d" % i for i in range(N)]
R = PolynomialRing(G, variables, order='lex')
R.inject_variables()

z = vector([eval(variables[i]) for i in range(N)])
t = randvect(G, N)
x = randvect(G, N)

f = z.dot_product(x)
ft = f(*t)
p = f - ft

dcmp = poly_decomp(p, z, t)
test_p = 0
for i, term in enumerate(dcmp):
    test_p += (eval(variables[i]) - t[i]) * term
assert(test_p == p)

print("Passed")
"""
