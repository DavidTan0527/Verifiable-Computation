import secrets

N = 5
UPPER_BOUND = 100
Fv = ZZ
x = vector(Fv, [secrets.randbelow(UPPER_BOUND) for _ in range(N)])
print(f"Verifiable Dot Product for dimension {N} vectors with {x}")
vdp = VerifiableDotProduct(x)

print("1. KeyGen with T = 32")
pk, sk = vdp.keygen(32)
print("G1 group size:", vdp.G1.order())
print("PK:", pk)
print("SK:", sk)
print()

# Define our dot product
def f(z):
    return sum(xx * zz for xx, zz in zip(x, z))

print("2. Setup function (public) key FK with PK and dot product function f")
fk = vdp.setup(pk)
print("FK:", fk)
print()

assert(fk.curve() == pk[5][0].curve()) # fk is in same group as g

print("3. Encrypt an input z with PK and a chosen random r")
print("x:", x)
z = vector(Fv, map(lambda xx : xx.strip(), input("z: ").split(",")))
c, Wc = vdp.encrypt(pk, z)
print("c:", c)
print()

print("4. Compute the dot product over the encrypted input z (also outputs signature)")
V, sgm = vdp.compute(pk, c, Wc)
print("V:", V)
print("σ:", sgm)
print()

print("5. Decrypt the encrypted output to get the result of dot product")
print("x * z:", f(z))
v = vdp.decrypt(sk, V)
print("v:", v)
print()
assert(f(z) == v)

print("6. Verify the signature")
valid = vdp.verify(pk, fk, z, v, sgm)
print("is valid?", valid)
assert(valid)
