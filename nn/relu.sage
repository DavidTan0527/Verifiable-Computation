from sage.structure.element import is_Vector

class VerifiableReLU(VerifiableNNLayer):
    """
    ReLU for values ranging between [-a, a] with a > 0 can be approximated by polynomial x^2 + a * x = 0
    The linear part can be scaled by c, i.e. ReLU(x) = cx, x > 0
    """
    def __init__(self, a, c = 1, T = 32, precision = 10000):
        self.a = a * precision
        self.c = c
        self.scale = self.c / (2 * self.a)
        self.precision = precision

        # make polynomial be non-negative
        self.f_min = Integer(ceil(self.a^2 / 4))

        if 2^T < self.a:
            print("Security parameter T too small for the chosen interval and precision")
            T = ceil(log(self.a, 2))
            print("Bumping T to", T)

        R.<x> = PolynomialRing(ZZ)
        self.protocol = VerifiablePolynomial(x^2 + self.a * x + self.f_min)
        pk, sk = self.protocol.keygen(T)

        self.pk = pk
        self._sk = sk

        fk = self.protocol.setup(pk)
        self.fk = fk

    def forward(self, x):
        if is_Vector(x):
            res = map(self._forward, x.list())
            return zip(*res)

        else:
            return self._forward(x)

    def _forward(self, x):
        xx = round(x * self.precision)
        if abs(xx) > self.a:
            raise ValueError(f"Input {x} is larger than configured interval [-{self.a/self.precision}, {self.a/self.precision}]")
        C = self.protocol.encrypt(self.pk, xx)
        V, sgm = self.protocol.compute(self.pk, C)
        y = self.protocol.decrypt(self._sk, V)
        y = float(y - self.f_min) / self.precision
        return y * self.scale, sgm

    def verify(self, pk, fk, x, y, sgm):
        if is_Vector(x):
            if not is_Vector(y) or not is_Vector(sgm):
                raise ValueError("x, y, sgm must all be vectors or all be non-vectors")

            res = map(lambda xx, yy, ss: self._verify(pk, fk, xx, yy, ss),
                      zip(x.list(), y.list(), sgm.list()))
            return all(res)

        else:
            if is_Vector(y) or is_Vector(sgm):
                raise ValueError("x, y, sgm must all be vectors or all be non-vectors")
            return self._verify(pk, fk, x, y, sgm)

    def _verify(self, pk, fk, x, y, sgm):
        x = Integer(round(x * self.precision))
        y = Integer(round(y * self.precision / self.scale))
        return self.protocol.verify(pk, fk, x, y, sgm)

