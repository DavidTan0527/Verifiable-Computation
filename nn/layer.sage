class VerifiableNNLayer:
    def forward(self, x):
        raise NotImplementedError()

    def verify(self, pk, fk, x, y, sgm):
        raise NotImplementedError()

