from pprint import pprint

G = ZZ
A = random_matrix(G, 4, 3)
dense = VerifiableDense(A.ncols(), A.nrows(), T=20, mat_precision=10, precision=10)
dense.set(A, vector([0]*4))
relu = VerifiableReLU(10, T=20, precision=10)

x = vector([0.4, 1.5, 1.7])
s = VerifiableSequential([dense, relu])

y, sgm = s.evaluate(x)
print(y)
pprint(sgm)

