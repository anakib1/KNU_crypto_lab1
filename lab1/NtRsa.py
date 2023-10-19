import random
from typing import Tuple
from NumberTheory import NumberTheory

class NTRSA:
    PRIMES = [414507281407, 2428095424619, 4952019383323, 12055296811267, 17461204521323, 53982894593057, 3657500101, 16123689073,88362852307,175650481151,414507281407]

    def __init__(self):
        self.p, self.q = random.sample(NTRSA.PRIMES, 2)
        self.n = self.p * self.q
        self.phi = (self.p - 1) * (self.q - 1)
        self.open_exp = random.sample([17, 257, 65537], 1)[0]
        if NumberTheory.gcd(self.open_exp, self.phi) != 1: 
            raise Exception("Exponent is not coprime. Should never happen, verify your configuration.")
        self.secret_exp = pow(self.open_exp, -1, self.phi)
        self.used_private = False

    def get_public_key(self) -> Tuple[int,int]:
        return (self.open_exp, self.n)
    
    def get_private_key(self) -> Tuple[int,int]:
        if self.used_private == True:
            return None
        self.used_private = True
        return (self.secret_exp, self.n)
    
class NTRSAClient:
    def encode(message:int, open_key:Tuple[int,int]) -> int:
        exp, n = open_key
        return pow(message, exp, n)
    def decode(message:int, secret_key:Tuple[int,int]) -> int:
        exp, n = secret_key
        return pow(message, exp, n)