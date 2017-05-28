import itertools
from Crypto.PublicKey import RSA
from Crypto.Cipher import PKCS1_OAEP

def coron_trivariate(pol, X, Y, Z, l=2, debug=False):
    P.<x,y,z> = PolynomialRing(ZZ)
    pol = pol(x,y,z)

    # Handle case where pol(0,0,0) == 0
    xoffset = 0

    while pol(xoffset,0,0) == 0:
        xoffset += 1

    pol = pol(x+xoffset,y,z)

    # Handle case where gcd(pol(0,0,0),X*Y*Z) != 1
    while gcd(pol(0,0,0), X) != 1:
        X = next_prime(X, proof=False)

    while gcd(pol(0,0,0), Y) != 1:
        Y = next_prime(Y, proof=False)

    while gcd(pol(0,0,0), Z) != 1:
        Z = next_prime(Z, proof=False)

    pol = P(pol/gcd(pol.coefficients())) # seems to be helpful
    p000 = pol(0,0,0)

    # maximum degree of any variable
    delta = max(pol.degree(x),pol.degree(y),pol.degree(z))

    W = max(abs(i) for i in pol(x*X,y*Y,z*Z).coefficients())
    u = W + ((1-W) % abs(p000))
    N = u*(X*Y*Z)^l # modulus for polynomials

    # Construct polynomials
    p000inv = inverse_mod(p000,N)
    polq = P(sum((i*p000inv % N)*j for i,j in zip(pol.coefficients(),
                                                 pol.monomials())))
    polynomials = []
    for i in range(delta+l+1):
        for j in range(delta+l+1):
            for k in range(delta+l+1):
                if 0 <= i <= l and 0 <= j <= l and 0 <= k <= l:
                    polynomials.append(polq * x^i * y^j * z^k * X^(l-i) * Y^(l-j) * Z^(l-k))
                else:
                    polynomials.append(x^i * y^j * z^k * N)

    # Make list of monomials for matrix indices
    monomials = []
    for i in polynomials:
        for j in i.monomials():
            if j not in monomials:
                monomials.append(j)
    monomials.sort()

    # Construct lattice spanned by polynomials with xX, yY, zZ
    L = matrix(ZZ,len(monomials))
    for i in range(len(monomials)):
        for j in range(len(monomials)):
            L[i,j] = polynomials[i](X*x,Y*y,Z*z).monomial_coefficient(monomials[j])

    # makes lattice upper triangular
    # probably not needed, but it makes debug output pretty
    L = matrix(ZZ,sorted(L,reverse=True))

    if debug:
        print "Bitlengths of matrix elements (before reduction):"
        print L.apply_map(lambda x: x.nbits()).str()
        set_verbose(2)

    L = L.LLL()

    if debug:
        print "Bitlengths of matrix elements (after reduction):"
        print L.apply_map(lambda x: x.nbits()).str()

    roots = []
    P2.<q> = PolynomialRing(ZZ)

    for i in range(L.nrows()-1):
        for j in range(i+1,L.nrows()):
            print "Trying rows %d, %d" % (i,j)
            pol2 = P(sum(map(mul, zip(L[i],monomials)))(x/X,y/Y,z/Z))
            pol3 = P(sum(map(mul, zip(L[j],monomials)))(x/X,y/Y,z/Z))

            r = pol.resultant(pol2, z)
            r2 = pol.resultant(pol3, z)
            r = r.resultant(r2,y)
            assert r.is_univariate()

            if r.is_constant(): # not independent
                continue

            r = r(q,0,0) # convert to univariate polynomial

            if len(r.roots()) > 0:
                for x0, _ in r.roots():
                    if x0 == 0:
                        continue
                    if debug:
                        print "Potential x0:",x0
                    for y0, _ in P2(r2(x0,q,0)).roots():
                        if debug:
                            print "Potential y0:",y0
                        for z0, _ in P2(pol(x0,y0,q)).roots():
                            if debug:
                                print "Potential z0:",z0
                            if pol(x0-xoffset,y0,z0) == 0:
                                roots += [(x0-xoffset,y0,z0)]
    return roots

pubkey = RSA.importKey(open("public.pem").read())
n,e = pubkey.n, pubkey.e
print(n)
print(e)
gamma = 0.4
delta = 0.1604

P.<x,y,z> = PolynomialRing(ZZ)
X = floor(n^delta)
Y = floor(n^(delta+1/2-gamma))
Z = floor(n^(delta+1/2-gamma))

pol = e^2*x^2 + e*x*(y+z-2) - (y+z-1) - (n-1)*y*z

roots = coron_trivariate(pol, X, Y, Z, l=0, debug=True)
if len(roots) > 0:
    print '!!!'
    x0,y0,z0 = roots[0]
    print(x0)
    print(y0)
    print(z0)
    d = int(x0)
    privkey = RSA.construct((n,e,d))
    cipher = PKCS1_OAEP.new(privkey)
    print(cipher.decrypt(open("flag.enc").read()))
