from pprint import pprint

G = ZZ  # work over the integers
T = 256 # security parameter
r = 4   # rank of the matrix
size = 6

A1 = random_matrix(G, r, size)
A2 = random_matrix(G, size, r)

model = VerifiableSequential([
    VerifiableDense(size, r, T=T, lbound=-1_000_000).set(A1),
    VerifiableDense(r, size, T=T, lbound=-1_000_000).set(A2),
    VerifiableDropout(0.1)
])

x = vector(random_matrix(G, size, 1).list())
y, sig = model.evaluate(x)

print(y)
print(A2 * A1 * x)

pk, fk = model.pk, model.fk

valid = model.verify(pk, fk, x, y, sig)
invalid = model.verify(pk, fk, x, vector(random_matrix(G, size, 1).list()), sig)
print(valid, invalid)
assert(valid and not invalid)

# G = ZZ
# dense = VerifiableDense(A.ncols(), A.nrows(), T=20, mat_precision=10, precision=10)
# dense.set(A, vector([0]*4))
# relu = VerifiableReLU(10, T=T, precision=10)

# x = vector([0.4, 1.5, 1.7])
# s = VerifiableSequential([dense, relu])

# y, sgm = s.evaluate(x)
# print(y)
# pprint(sgm)

