import multiprocess

def task(args):
    protocol = args[0]
    print(protocol)
    return True
    # return protocol.verify(*args)

class VerifiableMatMul:
    def __init__(self, A, b = None):
        if b is not None and A.base_ring() != b.base_ring():
            raise ValueError("A and b must have the same base ring")
        if b is not None and b.degree() != A.nrows():
            raise ValueError("A and b must have compatible sizes")

        self.A = A
        self.size_in, self.size_out = A.ncols(), A.nrows()
        self.b = b if b is not None else vector(A.base_ring(), [0] * self.size_out)
        self.protocols = [VerifiableDotProduct(A[i], self.b[i]) for i in range(self.size_out)]

    def keygen(self, T):
        return tuple(zip(*[protocol.keygen(T) for protocol in self.protocols]))

    def setup(self, pk):
        return tuple(protocol.setup(pki) for protocol, pki in zip(self.protocols, pk))

    def encrypt(self, pk, z, r):
        return tuple(protocol.encrypt(pki, z, r) for protocol, pki in zip(self.protocols, pk))

    def compute(self, pk, C):
        return tuple(zip(*[protocol.compute(pki, Ci) for protocol, pki, Ci in zip(self.protocols, pk, C)]))
        # with multiprocessing.Pool() as pool:
        #     args = zip(self.protocols, pk, C)
        #     _ = list(args)
        #     res = pool.starmap(task("compute"), args)
        #     return tuple(zip(*res))

    def decrypt(self, sk, V, bound):
        return tuple(protocol.decrypt(ski, Vi, bound) for protocol, ski, Vi in zip(self.protocols, sk, V))
        # with multiprocessing.Pool() as pool:
        #     args = zip(self.protocols, sk, V, [bound] * len(self.protocols))
        #     _ = list(args)
        #     res = pool.starmap(task("decrypt"), args)
        #     return tuple(res)

    def verify(self, pk, fk, z, v, sgm):
        return all(protocol.verify(pki, fki, z, vi, sgmi) for protocol, pki, fki, vi, sgmi in zip(self.protocols, pk, fk, v, sgm))
        # with multiprocess.Pool() as pool:
        #     args = zip(self.protocols, pk, fk, [z] * len(self.protocols), v, sgm)
        #     res = pool.map(task, args)
        #     # print(res)
        #     return all(res)


class MatMulWrapper:
    def __init__(self, A, b = None, T = 10, delta = 10_000, bound = 10_000_000):
        self.A = A
        self.b = vector(b) if b is not None else vector(A.base_ring(), [0] * A.nrows())

        self.delta = vector(b.base_ring(), [delta] * b.degree())
        self.bound = bound

        self.protocol = VerifiableMatMul(self.A, self.b + self.delta)
        pk, sk = self.protocol.keygen(T)

        self.pk = pk
        self._sk = sk

        fk = self.protocol.setup(pk)
        self.fk = fk

    def verify(self, x, y):
        return self.protocol.verify(self.pk, self.fk, x, y + self.delta, self.sgm)

    """
    Does A * x + b
    """
    def __call__(self, x, **kwargs):
        r = vector(random_matrix(x.base_ring(), 1, x.degree()).list())
        C = self.protocol.encrypt(self.pk, x, r)
        V, sgm = self.protocol.compute(self.pk, C)
        self.sgm = sgm
        return vector(self.protocol.decrypt(self._sk, V, self.bound)) - self.delta
