from timeit import timeit
"""
Ref: https://theory.stanford.edu/~dfreeman/papers/ants-embedding.pdf (Example 3.6)
EC with embedding degree 10
q = 6462310997348816962203124910505252082673338846966431201635262694402825461643
A = -3
B = 4946538166640251374274628820269694144249181776013154863288086212076808528141
n = 6462310997348816962203124910505252082512561846156628595562776459306292101261
k = 10

assert(q.is_prime())
F = Zmod(q).field() # Zmod ring promoted to field since q is prime

assert((q ^ k - 1) % n == 0)

Chose not to use because computation of torsion subgroup is extremely slow.
"""

"""
BN curves (https://eprint.iacr.org/2005/133.pdf, https://eprint.iacr.org/2010/429.pdf)
p = p(u) = 36u^4 + 36u^3 + 24u^2 + 6u + 1
n = n(u) = 36u^4 + 36u^3 + 18u^2 + 6u + 1
for some u ∈ Z.
"""
class Pairing:
    """
    Bilinear pairing using tate pairing on BN curves over GF
    """

    def __init__(self, uu=4965661367192848881):
        """
        uu: The parameter value to BN curve trace, curve order, and characteristic of Fp. This value can be any integer that gives a prime order.
        """
        # BN curve parameters
        t = 6*uu^2 + 1
        r = 36*uu^4 + 36*uu^3 + 18*uu^2 + 6*uu + 1
        p = r + t - 1
        assert(p.is_prime())
        assert(r.is_prime())
        assert(p % 4 == 3)
        k = 12

        self.p = p # base field prime order
        self.r = r # curve order
        self.k = k # embedding degree

        # https://hackmd.io/@vivi432/bn128-in-c
        Fp = GF(p)
        self.Fp = Fp

        P.<x> = Fp[]
        Fp2.<u> = Fp.extension(x^2 + 1)
        self.Fp2 = Fp2

        zeta = u+9
        assert(len(zeta.nth_root(2, all=True)) == 0
               and len(zeta.nth_root(3, all=True)) == 0)

        P.<t> = Fp2[]
        Fp12.<w> = Fp2.extension(t^6 - zeta)
        self.Fp12 = Fp12

        assert(Fp.order() ^ 2 == Fp2.order() and Fp2.order() ^ 6 == Fp12.order())

        b1 = 3
        self.G1 = EllipticCurve(Fp, [0, b1])
        assert(self.G1.order() == r)

        b2 = 3 * zeta^(-1)
        self.G2 = EllipticCurve(Fp2, [0, b2])
        assert(self.G2.order() % r == 0)

        self.G12 = EllipticCurve(Fp12, [0, b1])

        # Generators of the curves
        self.P = self.G1((1, 2))
        self.Q = self.G2((11559732032986387107991004021392285783925812861821192530917403151452391805634 * u
                          + 10857046999023057135944570762232829481370756359578518086990519993285655852781,
                           4082367875863433681332203403145435568316851327593401208105741076214120093531 * u
                          + 8495653923123431417604973247489272438418190587263600148770280649306958101930))

        self.m = r

        assert(self.m * self.P == 0 and self.m * self.Q == 0)

        # Used for point twisting
        self.w = Fp12([0, 1] + [0] * 10)

    # https://github.com/ethereum/py_ecc/blob/a1d18addb439d7659a9cbac861bf1518371f0afd/py_ecc/bn128/bn128_curve.py#L129
    def twist(self, P):
        _x, _y = P.xy()
        _x, _y = _x.polynomial().coefficients(), _y.polynomial().coefficients()

        _x = _x + [0] * (2 - len(_x))
        _y = _y + [0] * (2 - len(_y))

        xcoeffs = [_x[0] - _x[1] * 9, _x[1]]
        ycoeffs = [_y[0] - _y[1] * 9, _y[1]]

        nx = self.Fp12([xcoeffs[0]] + [0] * 5 + [xcoeffs[1]] + [0] * 5)
        ny = self.Fp12([ycoeffs[0]] + [0] * 5 + [ycoeffs[1]] + [0] * 5)

        return self.G12((nx * self.w^2, ny * self.w^3))

    def e(self, P, Q):
        if not (P.curve() == self.G1 and Q.curve() == self.G2):
            raise ValueError("Points do not lie on the curves defined")

        # Project both points to E(F_{p^12})
        Px = self.G12(P)
        Qx = self.twist(Q)

        assert(Px.parent() == Qx.parent())

        # TODO: issues parallelizing this
        return Px.tate_pairing(Qx, self.m, self.k, q=self.p)

    def test(self):
        m = self.m
        Px, Qx = self.P, self.Q

        print(f"Testing with basis ({Px}, {Qx})")

        print("--- order m ---")
        Px2_Qx3 = self.e(2*Px, 3*Qx)
        Px2_Qx3_m = Px2_Qx3 ^ m
        print("e(2*Px, 3*Qx):", Px2_Qx3)
        print("e(2*Px, 3*Qx)^m:", Px2_Qx3_m)
        assert(Px2_Qx3 != 1 and Px2_Qx3_m == 1)

        print("--- non-degeneracy ---")
        Px_Qx = self.e(Px, Qx)
        assert(Px_Qx != 1)
        print("e(Px, Qx) =", Px_Qx)
        N = 30
        print("Average time of pairing:", timeit("self.e(Px, Qx)", globals=dict(globals(), **locals()), number=N)/N)

        print("--- bilinearity ---")
        Px6_Qx7 = self.e(6*Px, 7*Qx)
        Px_Qx_42 = self.e(Px, Qx)^42
        assert(Px6_Qx7 == Px_Qx_42)
        print("Tested e(6*Px, 7*Qx) == e(Px, Qx)^42")

        Px6_Qx = self.e(6*Px, Qx)
        Px7_Qx = self.e(7*Px, Qx)

        assert(Px7_Qx == Px6_Qx * Px_Qx)
        print("Tested e(7*Px, Qx) == e(6*Px, Qx) * e(Px, Qx)")

        Px_Qx6 = self.e(Px, 6*Qx)
        Px_Qx7 = self.e(Px, 7*Qx)
        assert(Px_Qx7 == Px_Qx * Px_Qx6)
        print("Tested e(Px, 7*Qx) == e(Px, Qx) * e(Px, 6*Qx)")

