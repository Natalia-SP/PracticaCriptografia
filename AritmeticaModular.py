
# coding: utf-8

from random import randint
from functools import reduce

def ext_euclides(a, b):
    if b == 0:
        return a, 1, 0
    else:
        u1, u2, v1, v2 = 0, 1, 1, 0
        
        while b > 0: # Mientras que b sea mayor que 0
            q = a // b                                     # a = |_ a / b _|
            r, u, v = a - q * b, u2 - q * u1, v2 - q * v1  # r = a - qb, u = u2 - qu1, v = v2 -qv1
            a, b, u2, u1, v2, v1 = b, r, u1, u, v1, v      
            
        return a, u2, v2
    
def inverse(a,b):
    return ext_euclides(a,b)[1]


def big_pow(a, b, n):
    if b == 1:
        return a % n
    
    if a == 1:
        return 1
    
    a0, b0, p = a, b, 1
    
    while b0 > 0:
        # Si el bit está a 1, se incrementa el valor del producto
        if b0 % 2 == 1:
            p = p * a0
            
        b0, a0 = b0 // 2, a0*a0%n
            
    return p % n


def bifactor(num):
    a0, s = num, 0
    while a0 > 0 and a0 % 2 == 0:
        a0 //= 2
        s += 1
    return s, a0


def miller_rabin(num):
    # primos menores que 5
    if num == 2 or num == 3:
        return True
    # es 4 o menor que 4 o es par mayor que 2
    elif num == 4 or num < 2 or num % 2 == 0:
        return False
    else:
        s, u = bifactor(num - 1)   # Calculamos p-1 = 2^s * u
        a = randint(2, num - 2)    # a in [2, ..., p - 2]
        l = [big_pow(a, (2**i) * u, num) for i in range(s + 1)] # l = [a^u, a^2u, ..., a^2su] 
        # Primer elemento de la lista es 1 o -1
        if l[0] == 1 or l[0] == num - 1:
            return True # Probablemente primo

        # Ninguna de las potencias es igual a 1
        elif 1 not in l:
            return False # No es primo

        # Si aparece un 1 en la lista no precedido de un -1
        elif 1 in l and l[l.index(1)-1] != num - 1:
            return False # No primo

        # Si -1 está en la lista y no es el último elemento
        elif num - 1 in l and l[-1] != num - 1:
            return True  # Probablemente primo
    
    
def miller_rabin_test(num, n = 10):
    return reduce(lambda x, y: x and y, [miller_rabin(num) for i in range(n)])


def isqrt(n):
    x = n
    y = (x + 1) // 2
    while y < x:
        x = y
        y = (x + n // x) // 2
    return x

def baby_step_giant_step_original(a, b, p):
    # Si p es primo
    if miller_rabin_test(p):
        # Buscamos k tal que a^k = b, con a,b in Z_p
        if b == 1:
            return 0 # k = 0
        
        else:
            # Si k existe -> k = cs -r; 0 <= r < s; 1 <= c <= s
            s = isqrt(p - 1)
            # giant pass
            L = [pow(a, i*s, p) for i in range(1, s + 1)]
            # baby pass
            l = [(b * big_pow(a, i, p)) % p for i in range(s)]
            # calculamos la intersección entre L y l
            ks = list(filter(lambda x: x in L, l))
            # calculamos los k, que en caso de que p
            # no sea primitivo, habrá varios k
            for k in ks:
                yield (L.index(k) + 1) * s - l.index(k)
    
    else:
        print("p =", p, "no es primo.")
    
    return None


def baby_step_giant_step(a, b, p):
    # Si p es primo
    if miller_rabin_test(p):
        # Buscamos k tal que a^k = b, con a,b in Z_p
        if b == 1:
            return 0 # k = 0
        
        else:
            # Si k existe -> k = cs -r; 0 <= r < s; 1 <= c <= s
            # s = isqrt(p - 1)
            s = isqrt(p)
            # baby step
            l = {}
            for i in range(s):
                li = (b * big_pow(a, i, p)) % p
                l[li] = i
            # giants step
            for i in range(1, s + 1):
                # giant step
                Li = big_pow(a, i*s, p)
                # check if Li is on l
                li = l.get(Li)
                if li != None:
                    return i*s - li
            else:
                raise ValueError("No existe logaritmo para este número")
    else:
        raise AttributeError("p =", p, "no es primo.")


def Jacobi(a, p):
    if p % 2 != 0:
        symbol = 1  # inicializamos el símbolo de jacobi
        a0 = a % p  # 1: aplicamos (a / p) = (a % p / p)

        if a0 == 0:
            return 0
        elif a0 == 1:
            return 1
        elif a0 == -1:  # 5: -1 / p) = -1
            return ((-1) ** ((p - 1) // 2))

        u, s = bifactor(a0)  # 2: (ab / p) = (a / p)*(b / p)

        if u > 0:  # 3: (2 / p)  = (-1)**((p^2 - 1)/8)
            symbol = ((-1) ** ((p ** 2 - 1) // 8))
            if u % 2 == 0:
                symbol*=symbol


        # se puede descomponer n en a * b
        # y son distintos de 1 y -1
        # y p es impar
        if s == 1:  # 4: (1 / p)  = 1
            return symbol

        elif s == -1:  # 5: -1 / p) = -1
            return symbol * ((-1) ** ((p - 1) // 2))

        if p % 2 != 0:  # 6: (q / p)  = (-1)**((p - 1)(q - 1)/4) * (p / q)
            return symbol * Jacobi(p, s) * (-1) ** ((p - 1) * (s - 1) // 4)
        
    else:
        raise AttributeError('p tiene que ser impar')

        
def sqrt_mod(a, p):
    # si el número tiene raíz en Z_p
    if Jacobi(a, p) == 1:
        # buscamos un n tal que (n / p) = -1, es decir
        # n no es un residuo cuadrático
        for n in range(2, p):
            if Jacobi(n, p) == -1:
                break
        else:
            raise ValueError("Error: no tiene inverso")
            return None
        
        # Descomponemos p - 1 en 2^su
        u, s = bifactor(p - 1)
        # si u = 1, hacemos lo siguiente
        if u == 1:
            return big_pow(a, ((p+1) // 4), p)
        # en si u >= 2
        elif u >= 2:
            r, b, j, inv_a = big_pow(a, ((s+1) // 2), p), big_pow(n, s, p), 0, inverse(a, p)
            while j <= u - 2:
                if big_pow(inv_a*r**2, 2**(u - 2 - j), p) == p-1:
                    r *=b 
                    r %= p
                    
                b = b**2 
                j+=1
            
            return r
    else:
        raise ValueError("No tiene raíz cuadrada módulo", p)
        return None

# Devuelve las 2 raíces del número a en Z_p
def sqrts_mod(a, p):
    sqrt = sqrt_mod(a, p)
    sol = [sqrt, p - sqrt]
    sol.sort()
    return sol


def norm(cong):
    x, y, m = cong # separamos los coeficientes de la congruencia
    d, s = ext_euclides(x, m)[0:2] # obtenemos el mcd(x, m)
    # Si y mod d es distinto de 0, la congruencia no tiene solución
    if y % d != 0: 
        raise "Error: congruencia sin solución"
    # en caso contrario, normalizamos la congruencia
    else:
        # h = |_ y / d _|, f = |_ m / d _|
        h, f = y // d, m // d
        e = (h * s) % f
        return [1, e, f,]
    

def chinese_remainder(cong1, cong2, n):
    # cong1 => ax  = b  mod p => norm(cong1) => x = a mod p
    # cong2 => a'x = b' mod q => norm(cong2) => x = b mod q
    x1, a, p = norm(cong1)  # x1 = a + p*l
    x2, b, q = norm(cong2)  # a + p*l = b mod q
                            # p*l = (b - a) mod q
    inv_p = inverse(p, q)   # l = (b - a)*p^(-1) mod q
    aux = (b-a)*inv_p       # l = (b - a)*p^(-1) + q * s
    
    return (a + p*aux)%n    # x = a + p * ((b - a)*p^(-1))   ==> mod n


def sqrts_mod_n(a, p, q):
    sqrts_p, sqrts_q = sqrts_mod(a, p), sqrts_mod(a, q)
    n = p*q
    # Calculamos solo dos raíces, ya que las otras 
    # son las complementarias
    root1 = chinese_remainder([1, sqrts_p[0], p], [1, sqrts_q[0], q], n)
    root2 = chinese_remainder([1, sqrts_p[-1], p], [1, sqrts_q[0], q], n)
    return sorted([root1, n - root1, root2, n - root2])
    

def Fermat(n):
    #                       _       _
    # valor inicial de x = | sqrt(m) |
    #
    x = isqrt(n) + 1
    # mientras que val no sea un cuadrado perfecto
    while x < n:
        val = x**2 - n
        sqrt_val = isqrt(val)
        if val == sqrt_val**2:
            sqrt_val = int(sqrt_val)
            break
        x+=1
    else:
        # Si no lo encuentra, devuelve una lista vacía
        return []
    return [x - sqrt_val, x + sqrt_val]


def Pollard(n, c, f):
    # Inicializamos los elementos del algoritmo
    x, y, d = 2, 2, 1
    while d == 1:
        x, y = f(x), f(f(y))
        d = ext_euclides(abs(x - y), n)[0]
        print(x, "\t", y, "\t",d)
    
    
    if d == n:
        raise ValueError('No se encuentra descomposición para n con c =' + str(c))
    else:
        return d
    
    
def ρ_Pollard(n, c = 2):
    if c == 0 or c == -2:
        raise AttributeError('c tiene que ser distinto de 0 y -2')
    f = lambda x: (x**2 + c) % n
    factor = Pollard(n, c, f)
    return [factor, n // factor]

