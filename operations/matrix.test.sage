from sage.misc.profiler import Profiler

p = Profiler()

G = ZZ
A = random_matrix(G, 3, 3)
b = vector(random_matrix(G, A.nrows(), 1).list())
x = vector(random_matrix(G, A.ncols(), 1).list())

print(A); print()
print(x); print()
print(b); print()

print(A * x + b)

Av = VerifiableMatMul(A, b)
pk, sk = Av.keygen(32)
fk = Av.setup(pk)

C = Av.encrypt(pk, x)

p("calculate")
V, sgm = Av.compute(pk, C)
y = Av.decrypt(sk, V)

print("correct?", A * x + b == y)
p("verify")
print("verify?", Av.verify(pk, fk, x, y, sgm))

p()

print(p)

