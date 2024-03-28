
class VerifiableSequential:
    def __init__(self, layers : list[VerifiableNNLayer]):
        self._layers = layers


    def evaluate(self, x):
        sig = []
        y = x
        for layer in self._layers:
            y, sgm = layer.forward(y)
            sig.append(sgm)

        return y, sig

    def verify(self, pks, fks, x, y, sig):
        for layer, pk, fk, sgm in zip(self._layers, pks, fks, sig):
            inter_x, _ = layer.forward(x)
            if not layer.verify(pk, fk, x, inter_x, sgm):
                return False

            x = inter_x

        return x == y
