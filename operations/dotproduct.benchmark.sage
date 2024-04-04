import secrets
import time

from itertools import product

load("operations/dotproduct.paillier.sage")

Ns = [1, 2, 4, 8, 16, 32, 64]
Ts = [16, 32, 128]
rounds = 5
UPPER_BOUND = 100
Fv = ZZ

ntmp = {}
tnmp = {}

print(f"#rounds = {rounds}")
for N, T in product(Ns, Ts):
    total = 0
    for _ in range(rounds):
        x = vector(Fv, [secrets.randbelow(UPPER_BOUND) for _ in range(N)])
        vdp = VerifiableDotProduct(x)
        pk, sk = vdp.keygen(T)

        # Define our dot product
        def f(z):
            return sum(xx * zz for xx, zz in zip(x, z))

        fk = vdp.setup(pk)
        assert(fk.curve() == pk[5][0].curve()) # fk is in same group as g

        z = vector(Fv, [secrets.randbelow(UPPER_BOUND) for _ in range(N)])
        c = vdp.encrypt(pk, z)

        V, sgm = vdp.compute(pk, c)

        v = vdp.decrypt(sk, V)
        assert(f(z) == v)

        start = time.time()
        valid = vdp.verify(pk, fk, z, v, sgm)
        end = time.time()

        total += end - start

        assert(valid)

    print(f"N = {N}, T = {T:<3}, avg. time = {total/rounds:.5f}")
    if T not in tnmp: tnmp[T] = []
    if N not in ntmp: ntmp[N] = []
    tnmp[T].append((N, total/rounds))
    ntmp[N].append((T, total/rounds))

for T, arr in tnmp.items():
    print(T, "".join(map(str, arr)))
print("===")
for N, arr in ntmp.items():
    print(N, "".join(map(str, arr)))
