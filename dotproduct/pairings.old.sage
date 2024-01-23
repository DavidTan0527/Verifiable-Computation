from math import floor
from secrets import randbelow
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
p = p(u) = 36u^4 - 36u^3 + 24u^2 - 6u + 1
n = n(u) = 36u^4 - 36u^3 + 18u^2 - 6u + 1
for some u ∈ Z.
"""
class Pairing:
    """
    Bilinear pairing using weil pairing on BN curves over GF
    """

    def __init__(self):
        x = 4965661367192848881

        # BN cxrve parameters
        t = 6*x^2 + 1
        r = 36*x^4 + 36*x^3 + 18*x^2 + 6*x + 1
        p = r + t - 1
        assert(p.is_prime())
        assert(p % 4 == 3)
        k = 12
        self.k = k

        # https://hackmd.io/@vivi432/bn128-in-c
        Fp = GF(p)
        self.Fp = Fp

        P.<x> = Fp[]
        Fp2.<u> = Fp.extension(x^2 + 1)
        self.Fp2 = Fp2

        P.<y> = Fp2[]
        Fp12.<w> = Fp2.extension(y^2 - v)
        self.Fp12 = Fp12

        b1 = 3
        self.G1 = EllipticCurve(Fp, [0, b1])
        assert(self.G1.order() == r)

        b2 = 3*(9+u)^(-1)
        self.G2 = EllipticCurve(Fp2, [0, b2])
        assert(self.G2.order() % r == 0)

        # Generators of the curves
        self.P = self.G1((1, 2))
        self.Q = self.G2((11559732032986387107991004021392285783925812861821192530917403151452391805634 * u + 10857046999023057135944570762232829481370756359578518086990519993285655852781,
                           4082367875863433681332203403145435568316851327593401208105741076214120093531 * u + 8495653923123431417604973247489272438418190587263600148770280649306958101930))

        self.m = r

        assert(self.m * self.Px == 0 and self.m * self.Qx == 0)

        # try:
        #     self.Em = self.E
        # except ValueError:
        #     self.Em = self.E.change_ring(self.E.division_field(m))
        #     self.Px, self.Qx = self.Em.torsion_basis(m)

        # assert(self.Px.order() % m == 0 and self.Qx.order() % m == 0)
        # assert(m * self.Px == 0 and m * self.Qx == 0)

    m_q = 0xfffffffffffffffffffffffffffbffff
    m_a = ceil(log(m_q, 2))
    m_b = log(2^m_a - m_q - 1, 2)
    m_s = -1
    m_c = -1

    # B in E/F_p^2, A in E/F_p
    def EvalVertical1(self, B, A):
        (x_A, y_A) = A.xy()
        (x_B, y_B) = B.xy()
        r = x_B - x_A
        return r

    # B in E/F_p^2, A in E/F_p
    def EvalTangent1(self, B, A):
        (x_A, y_A) = A.xy()
        (x_B, y_B) = B.xy()
        if A == self.Fp(0):
            return 1
        if y_A == 0:
            return self.EvalVertical1(B, A)
        a = -3 * x_A^2
        b = 2 * y_A
        c = -b * y_A - a * x_A
        r = a * x_B + b * y_B + c
        return r

    # B in E/F_p^2, A1, A2 in E/F_p
    def EvalLine1(B, A1, A2):
        (x_A1, y_A1) = A1.xy()
        (x_A2, y_A2) = A2.xy()
        (x_B, y_B) = B.xy()
        if A1 == self.Fp(0):
            return self.EvalVertical1(B, A2)
        if A2 == self.Fp(0):
            return self.EvalVertical1(B, A1)
        if A1 == -A2:
            return self.EvalVertical1(B, A1)
        if A1 == A2:
            return self.EvalTangent1(B, A1)
        a = y_A1 - y_A2
        b = x_A2 - x_A1
        c = -b * y_A1 - a * x_A1
        r = a * x_B + b * y_B + c
        return r

    # B in E/F_p^2, A in E/F_p
    def TateMillerSolinas(A, B):
        v_num = self.Fp2(1)
        v_den = self.Fp2(1)
        t_num = self.Fp2(1)
        t_den = self.Fp2(1)
        V = A

        # Calculation of the (s * 2^b) contribution
        for i in range(b):
            t_num = t_num^2
            t_den = t_den^2
            t_num = t_num * self.EvalTangent1(B, V)
            V = 2*V
            t_den = t_den * self.EvalVertical1(B, V)

        # Normalization
        (x,y) = V.xy()
        V_b = self.Fp(x, s*y)

        # Accumulation
        if s == -1:
            v_num = v_num * t_den
            v_den = v_den * t_num * EvalVertical1(B, V)
        if s == 1:
            v_num = v_num * t_num
            v_den = v_den * t_den

        # Calculation of the 2^a contribution
        for i in range(b, a):
            t_num = t_num^2
            t_den = t_den^2
            t_num = t_num * EvalTangent1(B, V)
            V = 2*V
            t_den = t_den * EvalVertical1(B, V)

        # Normalization
        (x, y) = V.xy()
        V_a = EFp(x, s*y)

        # Accumulation
        v_num = v_num * t_num
        v_den = v_den * t_den

        # Correction for the (s * 2^b) and (c) contributions
        v_num = v_num * EvalLine1(B, V_a, V_b)
        v_den = v_den * EvalVertical1(B, V_a + V_b)

        if c == -1:
            v_den = v_den * EvalVertical1(B, A)

        # Correcting exponent
        eta = (p^2 - 1) / q

        r = (v_num / v_den)^eta
        return r
        
    # def twist(self, P):
    #     if P is None:
    #         return None
    #     _x, _y = P.xy()
    #     _xc, _yc = _x.coefficients(), _y.coefficients()

    #     xcoeffs = [_xc[0] - 9 * _xc[1], _xc[1]]
    #     ycoeffs = [_yc[0] - 9 * _yc[1], _yc[1]]

    #     nx = self.Fp12(xcoeffs[0] + xcoeffs[1] * w^6)
    #     ny = self.Fp12(ycoeffs[0] + ycoeffs[1] * w^6)


    def test(self):
        print(u, v, w)
        m = self.m
        Px, Qx = self.P, self.Q

        print(f"Testing with basis ({Px}, {Qx})")
        # assert(Px.order() == m and Qx.order() == m)

        print("--- order m ---")
        Px2_Qx3 = self.e(2*Px, 3*Qx)
        Px2_Qx3_m = Px2_Qx3 ^ m
        print("e(2*Px, 3*Qx):", Px2_Qx3)
        print("e(2*Px, 3*Qx)^m:", Px2_Qx3_m)
        assert(Px2_Qx3_m == 1)
        
        print("--- alternating ---")
        Px_Px = self.e(Px, Px)
        Qx_Qx = self.e(Qx, Qx)
        print("e(Px, Px):", self.e(Px, Px))
        print("e(Qx, Qx):", self.e(Qx, Qx))
        assert(Px_Px == 1)
        assert(Qx_Qx == 1)

        # print("--- non-degeneracy ---")
        Px_Qx = self.e(Px, Qx)
        # assert(Px_Qx != 1)
        # print("e(Px, Qx) =", Px_Qx)
        print("Average time of pairing:", timeit("self.e(Px, Qx)", globals=dict(globals(), **locals()), number=3000))

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



if __name__ == "__main__":
    # pr = Pairing(u=-(2^62 + 2^55 + 1), lm=10)
    pr = Pairing()
    pr.test()
