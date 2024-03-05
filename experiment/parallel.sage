import time
p = 103; A = 1; B = 18; E = EllipticCurve(GF(p), [A, B])
P = E(33, 91); n = P.order(); n
k = GF(n)(p).multiplicative_order(); k
P.tate_pairing(P, n, k)
Q = E(87, 51)

K.<a> = GF((p,k))
EK = E.base_extend(K); P = EK(P)
Qx = 69*a^5 + 96*a^4 + 22*a^3 + 86*a^2 + 6*a + 35
Qy = 34*a^5 + 24*a^4 + 16*a^3 + 41*a^2 + 4*a + 40
Q = EK(Qx, Qy);

h = 551269674; Q = h*Q
P = EK(P);

class Test:
    def pair_seq(self, P, Q, *args, **kwargs):
        return P.tate_pairing(Q, *args, **kwargs)

    @parallel(ncpus=16)
    def pair(self, P, Q, *args, **kwargs):
        return P.tate_pairing(Q, *args, **kwargs)


t = Test()
N = 100
Ps = [P] * N
Qs = [Q] * N
ns = [n] * N
ks = [k] * N

start = time.time()
res = list(t.pair(list(zip(Ps, Qs, ns, ks))))
print("Parallel:", time.time() - start)

start = time.time()
res = [t.pair_seq(P, Q, n, k) for P, Q, n, k in zip(Ps, Qs, ns, ks)]
print("Sequential:", time.time() - start)


