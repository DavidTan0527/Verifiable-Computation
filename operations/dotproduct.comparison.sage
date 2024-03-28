import secrets
import time

from itertools import product

Ts = [16, 32, 64, 128]
Ns = [1, 2, 4, 8, 16, 32]
bounds = [500, 2500, 5000]
# bounds = [500, 2500, 5000, 7500, 10000, 12500]
# UPPER_BOUND = 100
rounds = 5
Fv = ZZ

srcs = [
    "operations/dotproduct.sage",
    "operations/dotproduct.paillier.sage",
]

print(f"#rounds = {rounds}")
for N, T, upper_bound in product(Ns, Ts, bounds):
    totals = []
    for src in srcs:
        totals.append(0)
        for _ in range(rounds):
            load(src)
            x = vector(Fv, [secrets.randbelow(upper_bound) for _ in range(N)])
            vdp = VerifiableDotProduct(x)
            pk, sk = vdp.keygen(T)

            # Define our dot product
            def f(z):
                return sum(xx * zz for xx, zz in zip(x, z))

            fk = vdp.setup(pk)
            assert(fk.curve() == pk[5][0].curve()) # fk is in same group as g

            z = vector(Fv, [secrets.randbelow(upper_bound) for _ in range(N)])
            c = vdp.encrypt(pk, z)

            V, sgm = vdp.compute(pk, c)

            start = time.time()
            v = vdp.decrypt(sk, V, bound=N * upper_bound^2)
            end = time.time()
            assert(f(z) == v)

            # valid = vdp.verify(pk, fk, z, v, sgm)
            # assert(valid)

            totals[-1] += end - start


        
    total_b, total_p = totals
    print(f"N = {N}, T = {T}, Bound = {upper_bound}, BGN avg. = {total_b/rounds:<3.5f}, Paillier avg. = {total_p/rounds:<3.5f}")
