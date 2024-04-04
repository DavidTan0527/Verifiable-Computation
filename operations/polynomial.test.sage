from sage.misc.profiler import Profiler

p = Profiler()

q = 311308417372352325045323285119
assert(is_prime(q))

G = Zmod(q)
F.<x> = PolynomialRing(G)
f = 3 * x^2 + 203 * x + 5
vf = VerifiablePolynomial(f)

print("F:", F)
print("f:", f)

p("keygen")
pk, sk = vf.keygen(32)
print("pk:", pk)
print("sk:", sk)

p("setup")
fk = vf.setup(pk)
print("fk:", fk)

import secrets
z = secrets.randbelow(3000) - 1500

p("encrypt")
C = vf.encrypt(pk, z)
print("C:", C)

p("compute")
V, sgm = vf.compute(pk, C)
print("V:", V)
print("sgm:", sgm)

fz = f(z)
print("z", z)
print("f(z):", fz)

p("decrypt")
y = vf.decrypt(sk, V)
print("y", y)
print("Result correct?", y == fz)

p("verify")
verified = vf.verify(pk, fk, z, y, sgm)
print("Result verified?", verified)

p()

print(p)

