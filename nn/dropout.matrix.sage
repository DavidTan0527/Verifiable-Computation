import random

class VerifiableDropout(VerifiableNNLayer):
    def __init__(self, size, p = 0.5, T = 32, lbound = -10_000, precision = 10_000):
        self.size = size
        self.precision = precision
        self.delta = -lbound * precision if lbound < 0 else 0
        self.T = T
        
        if p < 0 or p > 1:
            raise ValueError(f"dropout probability has to be between 0 and 1, but got {p}")
        self.p = p

        self.protocol = VerifiableMatMul(identity_matrix(self.size),
                                         T = self.T,
                                         precision = self.precision,
                                         delta = self.delta)

        self.pk = self.protocol.pk
        self.fk = self.protocol.fk

    def _get_mask(self):
        mat = identity_matrix(self.size)
        for i in range(self.size):
            if random.random() < self.p:
                mat[i, i] = 0
        return mat

    def verify(self, pk, fk, x, y, sgm):
        return self.protocol.verify(pk, fk, x, y, sgm)

    def forward(self, x):
        if x.degree() != self.size:
            raise ValueError("dimension does not match")

        # Change to a different dropout mask
        mask = self._get_mask()
        self.protocol.updateFunc(A = mask)
        self.fk = self.protocol.fk

        y, sgm = self.protocol(x)
        return y, sgm


