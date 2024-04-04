from math import floor


class VerifiableDense(VerifiableMatMul, VerifiableNNLayer):
    def __init__(self, size_in, size_out, T = 32, lbound = -10_000, mat_precision = 1000, precision = 10_000):
        self.size_in, self.size_out = size_in, size_out
        self.mat_precision = mat_precision
        self.precision = precision
        self.T = T
        self.delta = -lbound * precision if lbound < 0 else 0

        # default to random weights and bias
        self.weights = random_matrix(RR, size_out, size_in)\
                            .apply_map(lambda x : floor(x * self.mat_precision))\
                            .change_ring(ZZ)

        self.bias = vector(random_matrix(RR, size_out, 1).list())\
                            .apply_map(lambda x : floor(x * self.mat_precision))\
                            .change_ring(ZZ)

        self.protocol = MatMulWrapper(self.weights,
                                      b = self.bias,
                                      T = self.T,
                                      precision = self.precision,
                                      delta = self.delta)

        self.pk = self.protocol.pk
        self.fk = self.protocol.fk

    def set(self, A = None, b = None):
        if A is not None:
            if A.ncols() != self.size_in or A.nrows() != self.size_out:
                raise ValueError("Dimension of A does not match")
            self.weights = A.apply_map(lambda x : floor(x * self.mat_precision))\
                            .change_ring(ZZ)
        if b is not None:
            if b.degree() != self.size_out:
                raise ValueError("Dimension of b does not match")
            self.bias = b

        self.protocol.updateFunc(self.weights, self.bias)

    def verify(self, pk, fk, x, y, sgm):
        y = y * self.mat_precision
        return self.protocol.verify(pk, fk, x, y, sgm)

    def forward(self, x):
        if x.degree() != self.size_in:
            raise ValueError("Dimension does not match")
        
        y, sgm = self.protocol(x)
        y = y / self.mat_precision
        return y, sgm


