# Warning: generation takes a while.

from Crypto.PublicKey import RSA
from Crypto.Cipher import PKCS1_OAEP

gamma = 0.4
delta = 0.1604
nbits = 4096

p = q = 0

g = random_prime(2^floor(nbits*gamma)-1, lbound=2^floor(nbits*gamma-1), proof=False)

while not p.is_prime(proof=False) or not q.is_prime(proof=False):
    count = 0

    while not p.is_prime(proof=False):
        a = randrange(2^floor(nbits*(1/2 - gamma)))
        p = 2*g*a + 1

    while not q.is_prime(proof=False):
        print "try",count
        b = randrange(2^floor(nbits*(1/2 - gamma)))
        while gcd(a,b) != 1 or not (2*g*a*b + a + b).is_prime(proof=False):
            b = randrange(2^floor(nbits*(1/2 - gamma)))
        q = 2*g*b + 1
        count += 1
        if count >= 25:
            break

assert p == 2*g*a + 1
assert q == 2*g*b + 1
assert p.is_prime(proof=False)
assert q.is_prime(proof=False)
assert g.is_prime(proof=False)
assert (2*g*a*b+a+b).is_prime(proof=False)
assert gcd(a,b) == 1

n = p*q
lcm = 2*g*a*b
d = e = k = 0

while gcd(k,2*g) != 1:
    d = randrange(2^floor(nbits*delta))
    while gcd(d,lcm) != 1:
        d = randrange(2^floor(nbits*delta))
    e = inverse_mod(d, lcm)
    k = (e*d-1)/lcm

n = long(n)
e = long(e)
d = long(d)
p = long(p)
q = long(q)

flag = open("flag.txt").read().strip()
key = RSA.construct((n,e,d,p,q))
cipher = PKCS1_OAEP.new(key)
ciphertext = cipher.encrypt(flag)
open("flag.enc","w").write(ciphertext)
open("private.pem","w").write(key.exportKey())
open("public.pem","w").write(key.publickey().exportKey())
