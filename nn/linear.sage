class VerifiableDense(VerifiableMatMul):
    def __init__(self, size_in, size_out, T = 32, lbound = -10_000, precision = 10_000):
        self.size_in, self.size_out = size_in, size_out
        self.precision = precision
        self.T = T
        self.delta = -lbound * precision if lbound < 0 else 0

        # default to random weights and bias
        self.weights = random_matrix(G, size_out, size_in)
        self.bias = vector(random_matrix(G, size_out, 1).list())

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
            self.weights = A
        if b is not None:
            if b.degree() != self.size_out:
                raise ValueError("Dimension of b does not match")
            self.bias = b

        self.protocol.updateFunc(self.weights, self.bias)

    def verify(self, pk, fk, x, y, sgm):
        return self.protocol.verify(pk, fk, x, y, sgm)

    def forward(self, x):
        if x.degree() != self.size_in:
            raise ValueError("Dimension does not match")
        
        y, sgm = self.protocol(x)
        return y, sgm


