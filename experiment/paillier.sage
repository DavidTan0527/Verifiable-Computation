T = 64
p = 0
q = 0

while p == q:
    p = random_prime(2^(T+1), False, lbound=2^T + 1)
    q = random_prime(2^(T+1), False, lbound=2^T + 1)

assert(is_prime(p) and is_prime(q) and p != q)

n = p * q
Zn = Zmod(n)
Zn2 = Zmod(n^2)
g = Zn2(n + 1)

assert(g.order() % n == 0)

def L(x):
    return (Zn2(x).lift() - 1) / n

lam = (p-1) * (q-1)
mu = 1 / Zn(lam)

pub = (n, g)
priv = (lam, mu)

print("public key:", pub)
print("private key:", priv)

# Encryption
def enc(m):
    r = randrange(1, n)
    while gcd(r, n) != 1:
        r = randrange(1, n)

    r = Zn2(r)
    c = g^m * r^n
    return c

m1 = randrange(0, n)
c1 = enc(m1)
print("m1:", m1)
print("c1:", c1)

m2 = randrange(0, n)
c2 = enc(m2)
print("m2:", m2)
print("c2:", c2)

c = c1 * c2
k = 5
ck = c1 ^ k
# Decryption
def dec(c):
    return Zn(L(c^lam)) * mu

res = dec(c)
print("decrypted m:", res)
print("correct?", res == m1 + m2)

print("const correct?", dec(ck) == k * m1)

