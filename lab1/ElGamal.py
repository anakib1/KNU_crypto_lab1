import collections
from NumberTheory import NTInverses
from typing import Tuple
from NumberTheory import IllegalArgumentException

def inv(n, q):
    return pow(n, -1, q)

Coord = collections.namedtuple("Coord", ["x", "y"])

# x^3 + ax + b = y^2
class EllipticCurve:
    def __init__(self, a:int, b:int, q:int):
        self.a = a
        self.b = b
        self.q = q
        self.zero = Coord(0, 0)
        pass

    def is_valid(self, p):
        if p == self.zero: return True
        l = (p.y ** 2) % self.q
        r = ((p.x ** 3) + self.a * p.x + self.b) % self.q
        return l == r

    def at(self, x : int):
        ysq = (x ** 3 + self.a * x + self.b) % self.q
        y, my = NTInverses.discrete_square_root(ysq, self.q)
        return Coord(x, y), Coord(x, my)

    def neg(self, p:Coord):
        return Coord(p.x, -p.y % self.q)

    def add(self, p1:Coord, p2:Coord):
        if p1 == self.zero: return p2
        if p2 == self.zero: return p1
        if p1.x == p2.x and (p1.y != p2.y or p1.y == 0):
            # p1 + -p1 == 0
            return self.zero
        if p1.x == p2.x:
            # p1 + p1: use tangent line of p1 as (p1,p1) line
            l = (3 * p1.x * p1.x + self.a) * inv(2 * p1.y, self.q) % self.q
            pass
        else:
            l = (p2.y - p1.y) * inv(p2.x - p1.x, self.q) % self.q
            pass
        x = (l * l - p1.x - p2.x) % self.q
        y = (l * (p1.x - x) - p1.y) % self.q
        return Coord(x, y)

    def mul(self, p, n):
        r = self.zero
        m2 = p
        # O(log2(n)) add
        while 0 < n:
            if n & 1 == 1:
                r = self.add(r, m2)
                pass
            n, m2 = n >> 1, self.add(m2, m2)
            pass
        return r

    def order(self, g):
        assert self.is_valid(g) and g != self.zero
        for i in range(1, self.q + 1):
            if self.mul(g, i) == self.zero:
                return i
        raise Exception("Invalid order")


class ElGamal(object):
    def __init__(self, curve : EllipticCurve, g:Coord, n = None):
        assert curve.is_valid(g)
        self.curve = curve
        self.g = g
        if n is not None:
            self.n = n
        else:
            self.n = curve.order(g)
        pass

    def generate_public_key(self, priv:int) -> Tuple[Coord, Coord]:
        return self.curve.mul(self.g, priv)

    def encode(self, plain:Coord, pub:Tuple[Coord, Coord], r : int) -> Tuple[Coord, Coord]:
        if not self.curve.is_valid(plain) or not self.curve.is_valid(pub):
            raise IllegalArgumentException("plain, pub must be points on a curve!")
        return (self.curve.mul(self.g, r), self.curve.add(plain, self.curve.mul(pub, r)))

    def decode(self, cipher:Tuple[Coord, Coord], priv:int) -> Coord:
        c1, c2 = cipher
        if not self.curve.is_valid(c1) and self.curve.is_valid(c2):
            raise IllegalArgumentException("c1, c2 must be points on a curve!!")
        return self.curve.add(c2, self.curve.neg(self.curve.mul(c1, priv)))

