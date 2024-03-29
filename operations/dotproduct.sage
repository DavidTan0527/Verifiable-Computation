from itertools import islice

def random_nonzero_element(G):
    return list(islice(filter(lambda x: x != 0 , (G.random_element() for _ in ZZ)), 1))[0]

class VerifiableDotProduct:
    def __init__(self, x, c = 0):
        """
        x: Dot product vector. Type: sage vector object
        c: Constant term added to the dot product. Type: base ring constant
        """
        self.dim = x.degree()
        self.x = x
        self.c = c

        # function setup
        variables = ["x%d" % i for i in range(self.dim)]
        R = PolynomialRing(self.x.base_ring(), variables, order='lex')
        R.inject_variables(verbose=False)

        z = vector([eval(variables[i]) for i in range(self.dim)])
        f = self.x.dot_product(z) + c
        self.f = f

    def keygen(self, T):
        """
        T: security parameter of the scheme

        Ran by model manager. Generates the key pair (PK, SK) with the parameters.
        """
        if T < 1:
            raise ValueError("Security parameter T must be a positive integer")

        pr = Pairing()
        G1 = pr.G1
        G2 = pr.G2
        G12 = pr.G12
        self.G1 = G1
        self.G2 = G2
        self.G12 = G12

        self.e = pr.e

        # Security requirement: q > 2^T
        q = random_prime(2^(T + 1), False, lbound=2^T + 1)
        self.q = q

        Zq = Zmod(q)
        s = random_nonzero_element(Zq)

        g = pr.P
        h = pr.Q

        gs = s * g
        g1 = G1.random_element()
        g2 = G1.random_element()

        # t = (t_1, t_2, ..., t_n) where t_i in Z*_q
        t = vector(Zq, [random_nonzero_element(Zq) for _ in range(self.dim)])

        W = [(t[i] * g1, t[i] * h) for i in range(self.dim)] + [sum(tt.lift()^2 for tt in t) * g1]

        PK = (q, G1, G2, G12, h, (g, gs), (g1, g2), W)
        SK = s

        # Note that in practice, the PK is not kept in the object, but passed in as it is public information.
        # We assign it to self to stick with the paper's function signatures
        self.pk = PK

        return PK, SK

    def setup(self, pk):
        """
        pk: Public Key generated by `keygen` method

        Ran by model manager. Setups the function `f` (dot product) with necessary algebraic structures.
        The function `f` is omitted from the parameters because it can be deduced from `x` in constructor.
        """
        q, G1, G2, G12, h, (g, gs), (g1, g2), W = pk
        assert(len(W) == len(self.x) + 1)

        fk = self.c * g1
        for (g1_t, h_t), xx in zip(W[:-1], self.x):
            fk += xx * g1_t

        # Add randomness
        # self.commit = vector(random_matrix(self.x.base_ring(), 1, self.dim))
        # fk += self.c * g1
        # for (g1_t, h_t), cc in zip(W[:-1], self.commit):
        #     fk += cc * g1_t

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

        q, G1, G2, G12, h, (g, gs), (g1, g2), W = pk

        r = vector(Zmod(q), (random_nonzero_element(Zmod(q)) for _ in range(self.dim)))

        # { g1^zi, g^ri, g2^zi, g^(s*ri) }
        C = [(z[i] * g1 + r[i] * g, z[i] * g2 + r[i] * gs) for i in range(self.dim)]

        return C

    def compute(self, pk, C):
        """
        pk: Public key generated by `keygen`
        C : Ciphertext of input to calculate dot product of

        Ran by the service provider. `f` omitted because it is implied in the `setup`
        """
        q, G1, G2, G12, h, (g, gs), (g1, g2), W = pk

        V = (self.c * g1, self.c * g2)

        for xx, cc in zip(self.x, C):
            V = (V[0] + xx * cc[0], V[1] + xx * cc[1])

        d = random_nonzero_element(Zmod(q))

        sgm0 = d * h
        self.sgm0 = sgm0

        sgm = []
        for (g1_t, h_t), xx in zip(W[:-1], self.x):
            sgm.append(d * (xx * g1 + g1_t))
        # for (g1_t, h_t), xx, cc in zip(W[:-1], self.x, self.commit):
        #     sgm.append(d * (xx * g1 + cc * g1 + g1_t))

        return (V, sgm)

    def decrypt(self, sk, V, bound=10_000_000):
        """
        sk: Secret Key generated by `keygen` method
        V : Encrypted output to be decrypted
        bound: Bound for solution to speed up discrete_log

        Ran by the model manager. Decrypts the computed (encrypted) output to the function output.
        """
        q, G1, G2, G12, h, (g, gs), (g1, g2), W = self.pk
        s = sk

        def pi(p): return s * p[0] - p[1]

        pi_V = pi(V)
        pi_g1_g2 = pi((g1, g2))

        # Find x such that x * pi((g1, g2)) == pi(V)
        return discrete_log_lambda(pi_V, pi_g1_g2, (0, bound), operation='+')

    def verify(self, pk, fk, z, v, sgm):
        q, G1, G2, G12, h, (g, gs), (g1, g2), W = pk

        assert(len(W) == len(z) + 1 and len(z) == len(sgm))

        H = W[-1]
        for (g1_t, h_t), zz in zip(W[:-1], z):
            H -= zz * g1_t

        verification = self.e(fk - v * g1 + H, self.sgm0)

        sgn = 1

        for (g1_t, h_t), zz, sgmi in zip(W[:-1], z, sgm):
            sgn *= self.e(sgmi, h_t - zz * h)

        # Check with pairing
        return verification == sgn

