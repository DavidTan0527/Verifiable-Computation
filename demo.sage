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
assert(valid and not invalid)

