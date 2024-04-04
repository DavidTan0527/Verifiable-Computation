
class VerifiableSequential:
    def __init__(self, layers : list[VerifiableNNLayer]):
        self._layers = layers
        self.fk = list(map(lambda p : p.fk, self._layers))
        self.pk = list(map(lambda p : p.pk, self._layers))


    def evaluate(self, x):
        sig = []
        y = x
        for layer in self._layers:
            y, sgm = layer.forward(y)
            sig.append((y, sgm))

        return y, list(zip(*sig))

    def verify(self, pks, fks, x, y, sig):
        inter_ys, sgms = sig
        if y != inter_ys[-1]:
            return False

        for layer, pk, fk, sgm, inter_y in zip(self._layers, pks, fks, sgms, inter_ys):
            print(layer)
            if not layer.verify(pk, fk, x, inter_y, sgm):
                return False
            x = inter_y

        return True
