class VerifiableReLU(VerifiablePolynomial):
    """
    ReLU for values ranging between [-a, a] with a > 0 can be approximated by polynomial x^2 + a * x = 0
    """

    def __init__(self, a = 1, F = ZZ):
        """
        a: Coefficient of linear term to fit the polynomial for different intervals. Type: base ring constant
        F: Field of the polynomial to work in.
        """
        R.<x> = PolynomialRing(F)
        super().__init__(x^2 + a * x)

class ReLUWrapper:
    def __init__(self, a, c = 1, T = 10, precision = 1000):
        self.a = a * precision
        self.c = c
        self.scale = self.c / (2 * self.a)
        self.precision = precision
        self.bound = ceil(2 * self.a / self.scale)

        if 2^T < self.bound:
            print("Security parameter T too small for the chosen interval and precision")
            T = ceil(log(self.bound, 2)) + ceil(log(self.a, 2))
            print("Bumping T to", T)

        self.protocol = VerifiableReLU(a = self.a)
        pk, sk = self.protocol.keygen(T)

        self.pk = pk
        self._sk = sk

        fk = self.protocol.setup(pk)
        self.fk = fk

    def forward(self, x):
        x = round(x * self.precision)
        C = self.protocol.encrypt(self.pk, x)
        V, sgm = self.protocol.compute(self.pk, C)
        y = self.protocol.decrypt(self._sk, V, self.bound)
        return float(y) * self.scale / self.precision, sgm

    def verify(self, pk, fk, x, y, sgm):
        x = round(x * self.precision)
        y = Integer(round(y * self.precision / self.scale))
        return self.protocol.verify(pk, fk, x, y, sgm)

