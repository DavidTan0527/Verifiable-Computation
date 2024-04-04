from math import floor


class VerifiableDense(VerifiableMatMul, VerifiableNNLayer):
    def __init__(self, size_in, size_out, T = 32, lbound = -10_000, mat_precision = 1000, precision = 10_000):
        self.size_in, self.size_out = size_in, size_out
        self.mat_precision = mat_precision
        self.precision = precision
        self.T = T
        self.delta = -lbound * precision if lbound < 0 else 0

        self.weights = random_matrix(ZZ, size_out, size_in)
        self.bias = vector(random_matrix(ZZ, size_out, 1).list())

        self.protocol = VerifiableMatMul(self.weights,
                                         b = self.bias,
                                         precision = self.precision,
                                         delta = self.delta)

        pk, sk = self.protocol.keygen(T)
        self.pk = pk
        self._sk = sk

        fk = self.protocol.setup(pk)
        self.fk = fk

    def set(self, A = None, b = None):
        if A is not None:
            if A.ncols() != self.size_in or A.nrows() != self.size_out:
                raise ValueError("Dimension of A does not match")
            self.weights = A

        if b is not None:
            if b.degree() != self.size_out:
                raise ValueError("Dimension of b does not match")
            self.bias = b
        else:
            self.bias = vector([0] * self.size_out)

        # Scale and cut off higher precision points
        self.weights = self.weights.apply_map(lambda x : floor(x * self.mat_precision))\
                            .change_ring(ZZ)

        self.bias = self.bias.apply_map(lambda x : floor(x * self.mat_precision))\
                             .change_ring(ZZ)

        self.protocol.updateFunc(self.weights, self.bias)
        fk = self.protocol.setup(self.pk)
        self.fk = fk

        return self

    def verify(self, pk, fk, x, y, sgm):
        y = y * self.mat_precision
        return self.protocol.verify(pk, fk, x, y, sgm)

    def forward(self, x):
        if x.degree() != self.size_in:
            raise ValueError("Dimension does not match")
        
        C = self.protocol.encrypt(self.pk, x)
        V, sgm = self.protocol.compute(self.pk, C)
        y = self.protocol.decrypt(self._sk, V)
        y = y / self.mat_precision
        return y, sgm


