import random

class VerifiableDropout(VerifiableNNLayer):
    def __init__(self, p = 0.5):
        if p < 0 or p > 1:
            raise ValueError(f"dropout probability has to be between 0 and 1, but got {p}")
        self.p = p

        self.pk = None
        self.fk = None

    def verify(self, pk, fk, x, y, sgm):
        num_zero = 0
        for xx, yy in zip(x, y):
            if yy == 0:
                num_zero += 1
            elif yy != xx:
                return False

        # p-value testing with some stats?
        return True

    def forward(self, x):
        y = [0 if random.random() < self.p else xx for xx in x]
        sgm = None
        return y, sgm


