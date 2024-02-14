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

    def encrypt(self, pk, z):
        return tuple(protocol.encrypt(pki, z) for protocol, pki in zip(self.protocols, pk))

    def compute(self, pk, C):
        return tuple(zip(*[protocol.compute(pki, Ci) for protocol, pki, Ci in zip(self.protocols, pk, C)]))

    def decrypt(self, sk, V, bound):
        return tuple(protocol.decrypt(ski, Vi, bound) for protocol, ski, Vi in zip(self.protocols, sk, V))

    def verify(self, pk, fk, z, v, sgm):
        return all(protocol.verify(pki, fki, z, vi, sgmi) for protocol, pki, fki, vi, sgmi in zip(self.protocols, pk, fk, v, sgm))
        # with multiprocess.Pool() as pool:
        #     # NOTE: works after removing self.protocols and pk (might be problem with pickling groups)
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

    def verify(self, pk, fk, x, y, sgm):
        return self.protocol.verify(pk, fk, x, y + self.delta, sgm)

    """
    Does A * x + b
    """
    def __call__(self, x, **kwargs):
        C = self.protocol.encrypt(self.pk, x)
        V, sgm = self.protocol.compute(self.pk, C)
        y = vector(self.protocol.decrypt(self._sk, V, self.bound)) - self.delta
        return y, sgm

