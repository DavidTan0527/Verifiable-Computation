G = ZZ
A = random_matrix(G, 32, 64)
# A = random_matrix(G, 8, 16)
b = vector(random_matrix(G, A.nrows(), 1).list())
x = vector(random_matrix(G, A.ncols(), 1).list())

print(A); print()
print(x); print()
print(b); print()

print(A * x + b)

Av = MatMulWrapper(A, b)
y = Av(x)

print(A * x + b == y)
print(Av.verify(x, y))

