T = 64
p = 0
q = 0

while p == q:
    p = random_prime(2^(T+1), False, lbound=2^T + 1)
    q = random_prime(2^(T+1), False, lbound=2^T + 1)

assert(is_prime(p) and is_prime(q) and p != q)

p = 17
q = 19

E = EllipticCurve(QQ, [1, -6])
Ep = E.change_ring(GF(p))
Eq = E.change_ring(GF(q))

EpEq = cartesian_product([Ep, Eq])

n = p * q
M = lcm(Ep.cardinality(), Eq.cardinality())
print(M)

Zn = Zmod(n)
Zn2 = Zmod(n^2)
E = E.change_ring(Zn2)
Q = E((1, 1, 52165)) # Find a way to select a random point?
# f = E.multiplication_by_m(n)

# print(f)
pub = (n, Q)
print(pub)
