class VerifiableMatMulInner:
    def __init__(self, A, b = None):
        if b is not None and b.degree() != A.nrows():
            raise ValueError("A and b must have compatible sizes")

        self.A = A
        self.size_in, self.size_out = A.ncols(), A.nrows()
        self.b = b if b is not None else vector(A.base_ring(), [0] * self.size_out)
        self.protocols = [VerifiableDotProduct(A[i], self.b[i]) for i in range(self.size_out)]

    def updateFunc(self, A, b = None):
        if self.A.dimensions() != A.dimensions() \
                or (b is not None and b.degree() != A.nrows()):
            raise ValueError("Invalid update to a different dimension")


        if b is not None:
            for i, (Ai, bi) in enumerate(zip(A.rows(), b)):
                self.protocols[i].updateFunc(Ai, bi)
        else:
            for i, Ai in enumerate(A.rows()):
                self.protocols[i].updateFunc(Ai)


    def keygen(self, T):
        return tuple(zip(*[protocol.keygen(T) for protocol in self.protocols]))

    def setup(self, pk):
        return tuple(protocol.setup(pki) for protocol, pki in zip(self.protocols, pk))

    def encrypt(self, pk, z):
        return tuple(protocol.encrypt(pki, z) for protocol, pki in zip(self.protocols, pk))

    def compute(self, pk, C):
        return tuple(zip(*[protocol.compute(pki, Ci) for protocol, pki, Ci in zip(self.protocols, pk, C)]))

    def decrypt(self, sk, V):
        return tuple(protocol.decrypt(ski, Vi).lift() for protocol, ski, Vi in zip(self.protocols, sk, V))

    def verify(self, pk, fk, z, v, sgm):
        return all(protocol.verify(pki, fki, z, vi, sgmi) for protocol, pki, fki, vi, sgmi in zip(self.protocols, pk, fk, v, sgm))


class VerifiableMatMul(VerifiableMatMulInner):
    def __init__(self, A, b = None, precision = 1, delta = 10_000):
        self.precision = precision

        self.A = A
        self.b = vector(b) if b is not None \
                    else vector(A.base_ring(), [0] * A.nrows())

        self.delta = vector(self.b.base_ring(), [delta] * self.b.degree())

        super().__init__(self.A, self.precision * (self.b + self.delta))

    def updateFunc(self, A = None, b = None):
        if A is not None:
            if A.nrows() != self.A.nrows() or A.ncols() != self.A.ncols():
                raise ValueError("dimension of A does not match")
            self.A = A

        if b is not None:
            if b.degree() != self.b.degree():
                raise ValueError("dimension of b does not match")
            self.b = b

        super().updateFunc(self.A, self.precision * (self.b + self.delta))

    def encrypt(self, pk, x):
        x = self.precision * x
        return super().encrypt(pk, x)

    def decrypt(self, sk, V):
        y = vector(super().decrypt(sk, V))
        y = y / self.precision - self.delta
        return y

    def verify(self, pk, fk, x, y, sgm):
        x = self.precision * x
        y = self.precision * (y + self.delta)
        return super().verify(pk, fk, x, y, sgm)
