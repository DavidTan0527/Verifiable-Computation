F.<x> = PolynomialRing(Zmod(10001))
f = 3 * x^2 + 203 * x + 10
vf = VerifiablePolynomial(f)

print("F:", F)
print("f:", f)

pk, sk = vf.keygen(12)
print("pk:", pk)
print("sk:", sk)

fk = vf.setup(pk)
print("fk:", fk)

z = 100
print(f"f({z}):", f(z))
C = vf.encrypt(pk, z)
print("C:", C)

V, sgm = vf.compute(pk, C)
print("V:", V)
print("sgm:", sgm)

y = vf.decrypt(sk, V, 1000000)
print("Result correct?", y == f(z))

verified = vf.verify(pk, fk, z, y, sgm)
print("Result verified?", verified)

