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
        return tuple(zip(*[protocol.compute(pki, Ci, Wci) for protocol, pki, (Ci, Wci) in zip(self.protocols, pk, C)]))

    def decrypt(self, sk, V):
        return tuple(protocol.decrypt(ski, Vi).lift() for protocol, ski, Vi in zip(self.protocols, sk, V))

    def verify(self, pk, fk, z, v, sgm):
        return all(protocol.verify(pki, fki, z, vi, sgmi) for protocol, pki, fki, vi, sgmi in zip(self.protocols, pk, fk, v, sgm))


class MatMulWrapper:
    def __init__(self, A, b = None, T = 10, precision = 1, delta = 10_000):
        self.precision_in = vector([precision] * A.ncols())
        self.inv_precision_in = vector([1/precision] * A.ncols())
        self.precision_out = vector([precision] * A.nrows())
        self.inv_precision_out = vector([1/precision] * A.nrows())

        self.A = A
        self.b = vector(b).pairwise_product(self.precision_out) if b is not None \
                else vector(A.base_ring(), [0] * A.nrows())

        self.delta = vector(self.b.base_ring(), [delta * precision] * self.b.degree())

        self.protocol = VerifiableMatMul(self.A, self.b + self.delta)
        pk, sk = self.protocol.keygen(T)

        self.pk = pk
        self._sk = sk

        fk = self.protocol.setup(pk)
        self.fk = fk

    def verify(self, pk, fk, x, y, sgm):
        x = x.pairwise_product(self.precision_in)
        y = y.pairwise_product(self.precision_out) + self.delta
        return self.protocol.verify(pk, fk, x, y, sgm)

    """
    Does A * x + b
    """
    def __call__(self, x, **kwargs):
        x = x.pairwise_product(self.precision_in)
        C = self.protocol.encrypt(self.pk, x)
        V, sgm = self.protocol.compute(self.pk, C)
        y = vector(self.protocol.decrypt(self._sk, V))
        y = (y - self.delta).pairwise_product(self.inv_precision_out)
        return y, sgm

