from math import ceil, floor, gcd
from itertools import groupby
from collections import Counter
from numpy import bitwise_xor, nonzero, bitwise_not, bitwise_and

def first_rule(seq):
    # Contamos el número de unos
    n_ones, n_bits = sum(seq), len(seq)
    # Si la longitud de la secuencia es 0
    if n_bits % 2 == 0:
        # a la fuerza tienen que ser iguales
        # num_ceros y num_unos
        return n_bits // 2 == n_ones
    else:
        n2_bits = n_bits / 2
        return ceil(n2_bits) == n_ones or floor(n2_bits) == n_ones


# Esta función calcula el número de rachas
# de longitud k en la secuencia seq      
def second_rule(seq):
    while seq[0] == seq[-1]:
        seq.append(seq.pop(0))
        
    # obtenemos todas las rachas posibles que existen en la secuencia
    runs = [list(g) for k, g in groupby(seq)]
    # y contamos el número de rachas que hay para cada longitud
    count = Counter(map(lambda x: len(x), runs))
    # comprobamos si se cumple que #runs(k+1) == runs(k)
    for i in range(1, len(count)):
        if count[i] != count[i+1]:
            if not count[i] >= 2*count[i+1]: # en el caso de que el siguiente elemento
                return False                 # no sea 1/2 veces más pequeño
    else:
        # si todo va bien, devolvemos True
        return True


def rotate(seq, length):
    seq = [seq[-1]] + seq[:length-1]


def third_rule(seq):
    length, rotated_seq = len(seq), seq
    rotate(rotated_seq, length)
    # Calculamos la distancia hamming como un xor entre
    # ambas cadenas, y contamos los bits no nulos
    norm = len(bitwise_xor(seq, rotated_seq).nonzero()[0])
    for i in range(1, length):
        rotate(rotated_seq, length)
        if norm != len(bitwise_xor(seq, rotated_seq).nonzero()[0]):
            return False
    else:
        return True


def Golomb(seq):
    if all(seq):
        return False
    else:
        rules = [first_rule, second_rule, third_rule]
        return all(map(lambda x: rules[x](seq), range(3)))
    

# Linear Feedback Shift Register
def LFSR(conex_poly, seed, length):
    register = list(seed)
    out = [0]*length
    # obtenemos el grado del polinomio de conexión
    for i in range(length):
        register.append((len(bitwise_and(conex_poly, register).nonzero()[0])) % 2)
        out[i] = register.pop(0)
    
    return out


f = lambda seed: ((seed[0] & seed[1]) | (not seed[2])) ^ seed[3]


def NLFSR(function, seed, length):
    out = [0]*length
    # obtenemos el grado del polinomio de conexión
    for i in range(length):
        seed.append(function(seed)) # añadimos el resultado al final
        out[i] = seed.pop(0) # y devolvemos el primer bit
    
    return out

# Convert a str to a binary's list like [0,1,0,0,1...]
def str_to_binlist(message):
    # get a binary list
    binary = list(bin(int.from_bytes(message.encode('utf-8'), byteorder='big'))) 
    # return the stream of 0 and 1
    binary.remove('b')
    return list(map(lambda x: int(x), binary))


def bin_to_str(bin_message):
    # pass elements from int to str 
    message = list(map(lambda x: str(x), bin_message))
    # create a full str
    complete_msg = ''.join(message)
    # convert to str using a integer in Z_2
    int_msg = int(complete_msg,2)
    return int_msg.to_bytes((int_msg.bit_length() + 7) // 8, 
                            byteorder='big').decode('utf-8', 'ignore')
    
    
def Geffe(arg1, arg2, arg3, length):
    # Calculate the LFSRs with their polynoms
    # and seeds with a fixed length
    LFSR1 = LFSR(arg1[0], arg1[1], length)
    LFSR2 = LFSR(arg2[0], arg2[1], length)
    LFSR3 = LFSR(arg3[0], arg3[1], length)
    # LFSR1 and LFSR2
    result1 = bitwise_and(LFSR1, LFSR2)
    # not LFSR2 and LFSR3
    result2 = bitwise_and(LFSR2, LFSR3)
    # result1 and result2
    return list(bitwise_xor(bitwise_xor(result1, result2), LFSR3))


# encrypt message using Geffe
def Geffe_encrypt(key, message):
    # get the binary list
    bin_message = str_to_binlist(message)
    # compute the Geffe sequence
    element = Geffe(key[0], key[1], key[2], len(bin_message))
    # encrypt the message
    cypher = bitwise_xor(element, bin_message)
    # return the cypher text
    return cypher


# decrypt message using Geffe 
def Geffe_decrypt(key, cypher):
    # compute the Geffe sequence
    element = Geffe(key[0], key[1], key[2], len(cypher))
    # Decrypt the message
    return bin_to_str(bitwise_xor(element,cypher))


def Berlekam_Massey(sequence):
    # initialize the parameters 
    n_bits = len(sequence)
    C, B, L, m = [0]*n_bits, [0]*n_bits, 0, -1
    C[0], B[0] = 1,1
    
    for N in range(n_bits):
        # Compute the next discrepancy
        di = (sequence[N] + sum([sequence[N - i]*C[ i] for i in range(1, L + 1)])) % 2
        if di == 1: # Update register
            T, D = C, list(B)
            # D = B * X^(N-m) (add zeros to left side)
            for i in range(N-m):
                D.insert(0,0); D.pop(-1) # remove the trash bits
            # C = C + B * X^(N-m)
            C = bitwise_xor(C, D)
            # Update the parameters
            if L <= N // 2:
                L, m, B = N +1 - L, N, T
                
    # return linear complexity and the polynom ready to use
    finalC = list(C[::-1])
    while finalC[0] == 0:
        finalC.pop(0)
    
    return L, finalC


lcm = lambda x, y: (x*y)//gcd(x, y)
