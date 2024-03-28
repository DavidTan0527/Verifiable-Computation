import time

from itertools import islice

def random_nonzero_element(G):
    return list(islice(filter(lambda x: x != 0 , (G.random_element() for _ in ZZ)), 1))[0]

class VerifiableDotProduct(SCC):
    def __init__(self, x, c = 0):
        """
        x: Dot product vector. Type: sage vector object
        c: Constant term added to the dot product. Type: base ring constant
        """
        self.dim = x.degree()
        self.updateFunc(x, c)

    def updateFunc(self, x, c = 0):
        if x.degree() != self.dim:
            raise ValueError("Invalid update to a different dimension")

        self.x = x
        self.c = c

        variables = ["x%d" % i for i in range(self.dim)]
        R = PolynomialRing(self.x.base_ring(), variables, order='lex')
        R.inject_variables(verbose=False)

        z = vector([eval(variables[i]) for i in range(self.dim)])

        f = self.x.dot_product(z) + c
        self.f = f

    def getW(self, g1, g, h):
        # t = (t_1, t_2, ..., t_n) where t_i in Z*_q
        t = vector(random_nonzero_element(Zmod(self.q)) for _ in range(self.dim))
        return [(t[i] * g1, t[i] * h) for i in range(self.dim)] + [sum(tt.lift()^2 for tt in t) * g1]

    # def keygen(self, T):
        """
        T: security parameter of the scheme

        Ran by model manager. Generates the key pair (PK, SK) with the parameters.
        """


    def setup(self, pk):
        """
        pk: Public Key generated by `keygen` method

        Ran by model manager. Setups the function `f` (dot product) with necessary algebraic structures.
        The function `f` is omitted from the parameters because it can be deduced from `x` in constructor.
        """
        q, G1, G2, G12, h, (g, gs), g1, W, (n, g_n2) = pk
        assert(len(W) == len(self.x) + 1)

        fk = self.c * g1
        for (g1_t, h_t), xx in zip(W[:-1], self.x):
            fk += xx * g1_t

        return fk

    def encrypt(self, pk, z):
        """
        pk: Public Key generated by `keygen` method
        z : Input to `f` (dot product). Type: sage vector object

        Ran by the client (independent execution from model manager and service provider).
        Encrypts `z` with `pk`.
        """
        if z.degree() != self.dim:
            raise ArithmeticError("Cannot encrypt vector of size %d (expected %d)" % (z.degree(), self.dim))

        if z.degree() != self.dim:
            raise ArithmeticError("Cannot use vector of size %d in current encryption setup (expected %d)" % (len(z), self.dim))

        return super().encrypt(pk, z)

    def compute(self, pk, C):
        """
        pk: Public key generated by `keygen`
        C : Ciphertext of input to calculate dot product of

        Ran by the service provider. `f` omitted because it is implied in the `setup`
        """
        q, G1, G2, G12, h, (g, gs), g1, W, (n, g_n2) = pk
        C, Wc = C

        r = random_nonzero_element(Zmod(n^2))
        V = g_n2 ^ self.c * r ^ n

        for xx, cc in zip(self.x, C):
            V *= cc ^ xx

        d = random_nonzero_element(Zmod(q))

        sgm0 = d * h
        self.sgm0 = sgm0

        sgm = []
        for (g1_t, h_t), xx in zip(W[:-1], self.x):
            sgm.append(d * (xx * g1 + g1_t))

        return (V, sgm)

    # def decrypt(self, sk, V):
        """
        sk: Secret Key generated by `keygen` method
        V : Encrypted output to be decrypted
        bound: Bound for solution to speed up discrete_log

        Ran by the model manager. Decrypts the computed (encrypted) output to the function output.
        """

    def verify(self, pk, fk, z, v, sgm):
        q, G1, G2, G12, h, (g, gs), g1, W, (n, g_n2) = pk

        assert(len(W) == len(z) + 1 and len(z) == len(sgm))

        H = W[-1]
        for (g1_t, h_t), zz in zip(W[:-1], z):
            H -= zz * g1_t

        verification = self.e(fk - v * g1 + H, self.sgm0)

        sgn = 1

        result = [self.e(sgmi, g1_t_h_t[1] - zz * h) for g1_t_h_t, zz, sgmi in zip(W[:-1], z, sgm)]
        for term in result:
            sgn *= term

        # Check with pairing
        return verification == sgn

