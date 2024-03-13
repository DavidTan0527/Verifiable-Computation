from itertools import islice

def random_nonzero_element(G):
    return list(islice(filter(lambda x: x != 0 , (G.random_element() for _ in ZZ)), 1))[0]

class MVP:
    def __init__(self, *args, **kwargs):
        raise NotImplementedError("Child classes must implement their own __init__")

    def getW(self, g1, h):
        raise NotImplementedError("Child classes must specify array W of supporting values for `keygen`")

    def keygen(self, T):
        """
        T: security parameter of the scheme

        Ran by model manager. Generates the key pair (PK, SK) with the parameters.
        """
        if T < 1:
            raise ValueError("Security parameter T must be a positive integer")

        pr = Pairing()
        self.pr = pr

        G1 = pr.G1
        G2 = pr.G2
        G12 = pr.G12
        self.G1 = G1
        self.G2 = G2
        self.G12 = G12

        self.e = pr.e

        q = random_prime(2^(T + 1), False, lbound=2^T + 1)
        self.q = q

        Zq = Zmod(q)
        s = random_nonzero_element(Zq)

        g = pr.P
        h = pr.Q

        gs = s * g
        g1 = G1.random_element()
        g2 = G1.random_element()

        W = self.getW(g1, g, h)

        # Paillier setup
        np = 0
        nq = 0

        while np == nq:
            np = random_prime(2^(T+1), False, lbound=2^T + 1)
            nq = random_prime(2^(T+1), False, lbound=2^T + 1)

        n = np * nq
        g_n2 = Zmod(n^2)(n + 1)

        lmd = (np-1) * (nq-1)
        mu = 1 / Zmod(n)(lmd)

        PK = (q, G1, G2, G12, h, (g, gs), (g1, g2), W, (n, g_n2))
        SK = (s, (lmd, mu))

        # Note that in practice, the PK is not kept in the object, but passed in as it is public information.
        # We assign it to self to stick with the paper's function signatures
        self.pk = PK

        return PK, SK

    def setup(self, pk):
        raise NotImplementedError("Child classes must implement their own `setup`")

    def encrypt(self, pk, z):
        q, G1, G2, G12, h, (g, gs), (g1, g2), W, (n, g_n2) = pk


        r = vector(random_nonzero_element(Zmod(n^2)) for _ in range(len(z)))
        # Client is able to retrieve the random vector for verification later
        self.r = r

        C = [g_n2 ^ z[i] * r[i] ^ n for i in range(len(z))]
        Wc = [z[i] * g1 + r[i] * g for i in range(len(z))]

        return C, Wc

    def compute(self, pk, C):
        raise NotImplementedError("Child classes must implement their own `compute`")

    def decrypt(self, sk, V):
        q, G1, G2, G12, h, (g, gs), (g1, g2), W, (n, g_n2) = self.pk
        s, (lmd, mu) = sk

        def L(x): return (Zmod(n^2)(x).lift() - 1) / n
        
        return Zmod(n)(L(V^lmd)) * mu

    def verify(self, pk, fk, z, v, sgm):
        raise NotImplementedError("Child classes must implement their own `verify`")

