from sage.misc.profiler import Profiler

p = Profiler()

G = ZZ
A = random_matrix(G, 16, 32)
# A = random_matrix(G, 2, 3)
b = vector(random_matrix(G, A.nrows(), 1).list())
x = vector(random_matrix(G, A.ncols(), 1).list())

print(A); print()
print(x); print()
print(b); print()

print(A * x + b)

Av = MatMulWrapper(A, b)

p("calculate")
y, sgm = Av(x)

pk, fk = Av.pk, Av.fk

print("A * x + b =", y)

print("correct?", A * x + b == y)
p("verify")
print("verify?", Av.verify(pk, fk, x, y, sgm))

p()

print(p)

