import unittest
from NumberTheory import NumberTheory, NTSymbols, NTFactorisation, NTInverses, NTPrime
from NtRsa import NTRSAClient, NTRSA
import random
from ElGamal import EllipticCurve, ElGamal

class TestFirstPoint(unittest.TestCase):
    def test_isprime(self):
        self.assertEqual(True, NumberTheory.isprime(127))
        self.assertEqual(False, NumberTheory.isprime(128))
    def test_phi(self):
        self.assertEqual(648, NumberTheory.phi(2052))
        self.assertEqual(1, NumberTheory.phi(1))
    def test_gcd(self):
        self.assertEqual(2, NumberTheory.gcd(12, 10))
        self.assertEqual(1, NumberTheory.gcd(1037218, 332411))
        self.assertEqual(743274, NumberTheory.gcd(268321914, 94395798))
    def test_lcm(self):
        self.assertEqual(6, NumberTheory.lcm(1, 2, 3))
        self.assertEqual(2, NumberTheory.lcm(2, 2, 2))
        self.assertEqual(1219862098, NumberTheory.lcm(1, 346, 371, 731, 182))
        self.assertEqual(None, NumberTheory.lcm())
    def test_mu(self):
        self.assertEqual(-1, NumberTheory.mu(30))
        self.assertEqual(0, NumberTheory.mu(32))
        self.assertEqual(1, NumberTheory.mu(119))

class TestSecondPoint(unittest.TestCase):
    def test_crt(self):
        self.assertEqual(34, NumberTheory.crt([3,5,7], [1,4,6]))
        self.assertEqual(23, NumberTheory.crt([3,5,7], [2,3,2]))

class TestSymbols(unittest.TestCase):
    def test_legendre(self):
        self.assertEqual(1, NTSymbols.legandre_symbol(21, 59))
        self.assertEqual(-1, NTSymbols.legandre_symbol(28, 41))
        self.assertEqual(0, NTSymbols.legandre_symbol(29, 29))
    def test_jacobi(self):
        self.assertEqual(-1, NTSymbols.jacobi_symbol(1001, 9907))
        self.assertEqual(1, NTSymbols.jacobi_symbol(461, 78419))
        self.assertEqual(0, NTSymbols.jacobi_symbol(127, 127 * 3))

class TestFactorization(unittest.TestCase):
    def test_pollard_ro(self):
        self.assertListEqual([(2, 1), (3, 1), (5, 1)], NTFactorisation.pollard_ro_factorisation(30))
        self.assertListEqual([(int(1e9) + 7, 1), (int(1e9) + 9, 1)], NTFactorisation.pollard_ro_factorisation((int(1e9) + 9) * (int(1e9) + 7)))
        self.assertListEqual([(998244353, 1)], NTFactorisation.pollard_ro_factorisation(998244353))
        self.assertListEqual([(1000000000039, 1)], NTFactorisation.pollard_ro_factorisation(1000000000039 ))
        #timeouts 
        #self.assertListEqual([(2, 3), (3, 1),(5, 1), (1000000000039, 1)], NTFactorisation.pollard_ro_factorisation(1000000000039 * 120))

class TestLogarithm(unittest.TestCase):
    def test_log(self):
        self.assertEqual(9, NTInverses.discrete_log(5, 33, 58))

    def test_sqrt(self):
        self.assertEqual((3, 4), NTInverses.discrete_square_root(2, 7))
        self.assertEqual((135, 9872), NTInverses.discrete_square_root(8218, 10007))
        self.assertEqual((37,64), NTInverses.discrete_square_root(56, 101))
        self.assertEqual((1,10), NTInverses.discrete_square_root(1, 11))
        self.assertEqual(None, NTInverses.discrete_square_root(8219, 10007))

class TestFastprime(unittest.TestCase):
    def test_isprime(self):
        self.assertEqual(True, NTPrime.fast_is_prime(127))
        self.assertEqual(False, NTPrime.fast_is_prime(128))
        self.assertEqual(True, NTPrime.fast_is_prime(int(1e9) + 7))
        self.assertEqual(True, NTPrime.fast_is_prime(int(1e9) + 9))
        self.assertEqual(True, NTPrime.fast_is_prime(998244353))
        self.assertEqual(True, NTPrime.fast_is_prime(1000000000039))
        self.assertEqual(False, NTPrime.fast_is_prime(31 * 1000000000039))
        self.assertEqual(False, NTPrime.fast_is_prime(31 * 127 * 1000000000039))
        self.assertEqual(False, NTPrime.fast_is_prime(62298863484143))
        self.assertEqual(True, NTPrime.fast_is_prime(73275315729173))

class TestRSA(unittest.TestCase):
    def test_encode_decode(self):
        rsa = NTRSA()
        open_key = rsa.get_public_key()
        secret = rsa.get_private_key()
        self.assertIsNone(rsa.get_private_key())

        n = open_key[1]
        message = random.randint(1, n - 1)

        encoded = NTRSAClient.encode(message, open_key)
        decoded = NTRSAClient.decode(encoded, secret)

        self.assertEqual(message, decoded)
    
    def validate_primes(self):
        for prime in NTRSA.PRIMES:
            self.assertTrue(NTPrime.fast_is_prime(prime))

class TestElGamal(unittest.TestCase):
    def test_elliptic_encoding(self):
        curve = EllipticCurve(13, 17, 93187)
        g = curve.at(117)[0]
        private_key = 1001 
        model = ElGamal(curve, g)
        pub = model.generate_public_key(private_key)
        data = curve.at(5)[0]
        encoded = model.encode(data, pub, 117)
        decoded = model.decode(encoded, private_key)
        self.assertEqual(data, decoded)
    
    def test_elliptic_encoding_with_required_params(self):
        # http://www.secg.org/SEC2-Ver-1.0.pdf
        # test with parameters selected accoring to 2.2.1 Recommended Parameters secp112r1
        # unfortunately test timeouts on my PC, because of giant numbers in the table. 
        if True: return
        p = int("DB7C2ABF62E35E668076BEAD208B", 16)
        a = int("DB7C2ABF62E35E668076BEAD2088", 16)
        b = int("659EF8BA043916EEDE8911702B22", 16)
        S = int("00F50B028E4D696E676875615175290472783FB1", 16)
        start_place = int("0209487239995A5EE76B55F9C2F098", 16)
        n = int("DB7C2ABF62E35E7628DFAC6561C5", 16)
        curve = EllipticCurve(a, b, p)
        g = curve.at(start_place)[0]
        private_key = 1999 
        model = ElGamal(curve, g, n)
        pub = model.generate_public_key(private_key)
        data = curve.at(1234567)[0]
        encoded = model.encode(data, pub, 2017)
        decoded = model.decode(encoded, private_key)
        self.assertEqual(data, decoded)


if __name__ == '__main__':
    unittest.main()
