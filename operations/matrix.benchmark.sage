import secrets
import time

from itertools import product

load("operations/matrix.sage")

Ms = [2, 4, 8, 16, 32, 48, 64]
Ns = [2, 4, 8, 16, 32, 48, 64]
rounds = 5
UPPER_BOUND = 100
G = ZZ

ntmp = {}
tnmp = {}

print(f"#rounds = {rounds}")
for M, N in list(product(Ms, Ns))[34:]:
    comp_total = 0
    verf_total = 0
    for _ in range(rounds):
        A = random_matrix(Zmod(10000), M, N).lift()
        b = vector(random_matrix(Zmod(10000), A.nrows(), 1).list()).lift()
        x = vector(random_matrix(Zmod(10000), A.ncols(), 1).list()).lift()

        Av = VerifiableMatMul(A, b, T = 128)
        pk, fk = Av.pk, Av.fk


        start = time.time()
        y, sgm = Av(x)
        end = time.time()
        if A * x + b != y:
            print("ERROR", "answer wrong")
            print("Expected:", A*x + b)
            print("Actual:", y)

        comp_total += end - start

        start = time.time()
        valid = Av.verify(pk, fk, x, y, sgm)
        end = time.time()
        if not valid:
            print("ERROR", "verification wrong")

        verf_total += end - start

        assert(valid)

    print(f"size: {M:2} x {N:2} , avg. (compute, verify) = {comp_total/rounds:.5f}, {verf_total/rounds:.5f}")
    # if T not in tnmp: tnmp[T] = []
    # if N not in ntmp: ntmp[N] = []
    # tnmp[T].append((N, total/rounds))
    # ntmp[N].append((T, total/rounds))

# for T, arr in tnmp.items():
    # print(T, "".join(map(str, arr)))
# print("===")
# for N, arr in ntmp.items():
    # print(N, "".join(map(str, arr)))
