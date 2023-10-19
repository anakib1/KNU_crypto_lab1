#region imports

from typing import List, Tuple
import math
from functools import reduce
from itertools import product
import random 

#endregion imports

# region utils

class IllegalArgumentException(Exception):
    pass

def integralFunction(func):
    def wrapper(*args, **kwargs):
        for arg in args:
            if not isinstance(arg, int):
                raise IllegalArgumentException(f"Argument '{arg}' is not of type int")
        return func(*args, **kwargs)
    return wrapper

# endregion utils

class NumberTheory:
    @integralFunction
    def isprime(n: int) -> bool:
        d = 2
        while d * d <= n:
            if n % d == 0:
                return False 
            d += 1
        return True

    @integralFunction
    def phi(n : int) -> int:
        ans = n
        d = 2
        while d * d <= n:
            if n % d == 0:
                while n % d == 0: 
                    n //= d
                ans -= ans // d
            d += 1
        if n > 1: ans -= ans // n
        return ans

    @integralFunction
    def gcd(x:int, y:int) -> int:
        while y != 0:
            x, y = y, x % y
        return x

    @integralFunction
    def lcm(*args:List[int]):
        if len(args) == 0:
            return None 
        x = args[0]
        for y in args[1:]:
            x = x // NumberTheory.gcd(x, y) * y 
        return x 
    
    @integralFunction
    def mu(n:int) -> int:
        p = 0
        if n % 2 == 0:
            n //= 2
            p += 1
            if n % 2 == 0:
                return 0
        
        d = 3
        while d * d <= n: 
            if n % d == 0:
                n /= d
                p += 1
                if n % d == 0:
                    return 0
            d += 2
        
        return (p % 2) * 2 - 1 
    
    def crt(m:List[int], r:List[int]) -> int:
        for i in range(len(m)):
            for j in range(i + 1, len(m)):
                if NumberTheory.gcd(m[i], m[j]) > 1:
                    raise NotImplementedError(f"CRT for not coprime modulos haven't bin implemented yet. ({m[i]}, {m[j]}) > 1")
        sum = 0 
        M = reduce(lambda x,y: x * y, m)
        for _m, _r in zip(m, r):
            p = M // _m 
            sum += _r * pow(p, -1, _m) * p
        return sum % M 
    
class NTSymbols:
    @integralFunction
    def legandre_symbol(a:int, p:int) -> int:
        a %= p
        if p == 2 or not NumberTheory.isprime(p):
            raise IllegalArgumentException("P must be prime, got: " + str(p))
        res = pow(a, (p - 1)//2, p)
        if res == p - 1:
            return -1
        return res
    
    @integralFunction
    def jacobi_symbol(n : int, k : int) -> int:
        if k % 2 != 1:
            raise IllegalArgumentException("K must be odd, got: " + str(k))
        n = n % k
        t = 1
        while n != 0:
            while n % 2 == 0:
                n //= 2
                r = k % 8
                if r == 3 or r == 5: t = -t
            n, k = k, n
            if n % 4 == 3 and k % 4 == 3: t = -t
            n = n % k
        return t if k == 1 else 0

class NTFactorisation:

    USUAL_BOUND = int(1e6)

    @integralFunction
    def pollard_ro_factorisation(n : int) -> List[Tuple[int, int]]:

        def g(x:int, m : int) -> int:
            return (x * x + 1) % m
        
        def run_ro(start_value : int, m : int) -> int:
            x = start_value
            y = start_value
            d = 1
            while d == 1:
                x = g(x, m)
                y = g(g(y, m), m)
                d = NumberTheory.gcd(abs(x - y), m)
            if d == m:
                return None 
            return d
        
        def naive(m : int ) -> List[Tuple[int, int]]:
            ret = []
            d = 2
            while d * d <= m:
                if m % d == 0:
                    cnt = 0
                    while m % d == 0:
                        cnt += 1
                        m /= d
                    ret.append((d, cnt))
                d += 1
            
            if m > 1:
                ret.append((m, 1))

            return ret

        factors = []
        to_factorise = [n]
        while len(to_factorise) > 0:
            m = to_factorise[0]
            if m > NTFactorisation.USUAL_BOUND:
                if NTPrime.fast_is_prime(m):
                    factors.append((m, 1))
                else:
                    ret = None
                    for value in [2, 3, 7, 19, 119, 127, 37, 34, 88]:
                        ret = run_ro(value, m)
                        if ret != None:
                            break
                    if ret == None:
                        # probaby prime but warn...
                        factors.append((m, 1))
                    else:
                        to_factorise += [ret, n//ret]
            else:
                factors += naive(m)
            to_factorise.pop(0)

        sums = {}
        for (x, y) in factors:
            sums[x] = sums.setdefault(x, 0) + y  

        return sorted([(x, sums[x]) for x in sums.keys()])

class NTInverses:
    """
    a^x === b mod p
    """
    @integralFunction
    def discrete_log(a : int, b : int, p : int) -> int:
        m = math.ceil(math.sqrt(p - 1))
        arr = {pow(a, j, p) : j  for j in range(m)}
        ainversem = pow(a, m * (NumberTheory.phi(p) - 1), p)
        y = b
        for j in range(m):
            if y in arr:
                return j * m + arr[y]
            y *= ainversem
            y %= p
        return None

    @integralFunction
    def discrete_square_root(n : int, p : int) -> int: 
        n %= p
        if n == 0 or n == 1:
            return (n, -n % p)
        if NTSymbols.legandre_symbol(n, p) != 1:
            return None
        
        if p % 4 == 3:
            ans = pow(n,(p+1)//4,p)
            return tuple(sorted([ans, -ans % p]))
        
        good = False
        for a in range(1, p):
            if NTSymbols.legandre_symbol(a * a - n, p) == -1:
                good = True
                break
        if not good: return None

        w = (a * a - n) % p
        
        x = (a, 1)

        def my_prod(u, v):
            return ((u[0] * v[0] + w * u[1] * v[1]) % p, (u[0] * v[1] + u[1] * v[0]) % p) 
        
        def fast_pw(x, y):
            if y == 0:
                return (1, 0)
            xx = fast_pw(x, y // 2)
            z = my_prod(xx, xx)
            if (y % 2 == 1):
                z = my_prod(z, x)
            return z
        
        res = fast_pw(x, (p+1)//2)
        return tuple(sorted([res[0], p - res[0]]))

class NTPrime:
    BASES = [2, 3, 5, 7, 11]
    def fast_is_prime(n : int) -> bool:
        if n <= 1: return False
        if n == 2: return True
        if n % 2 == 0 : return False

        d = n - 1
        s = 0
        while d % 2 == 0: 
            s += 1
            d //= 2
        for a in NTPrime.BASES:
            x = pow(a, d, n)
            for i in range(s):
                y = (x * x) % n
                if y == 1 and x != 1 and x != n - 1:
                    return False
                x = y
            if y != 1:
                return False
            
        return True
