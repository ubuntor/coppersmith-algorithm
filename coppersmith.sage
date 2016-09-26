def coron(pol, X, Y, k=2):
    """
    Applies Coron's reformulation of Coppersmith's algorithm for finding small
    integer roots of bivariate polynomials modulo an integer.

    pol: The polynomial to find small integer roots of.
    X: Upper limit on x.
    Y: Upper limit on y.
    k: Determines size of lattice. Increase if this fails.
    """

    if pol.nvariables() != 2:
        raise ValueError("pol is not bivariate")

    if pol(0,0) == 0:
        raise ValueError("pol(0,0) == 0 not supported (yet)")

    if gcd(pol(0,0), X*Y) != 1:
        raise ValueError("gcd(pol(0,0), X*Y) != 1 not supported (yet)")

    P.<x,y> = PolynomialRing(ZZ)
    pol = pol(x,y)
    pol = P(pol/gcd(pol.coefficients())) # seems to be helpful

    coeffs = [pol.monomial_coefficient(P(1)), pol.monomial_coefficient(x),
              pol.monomial_coefficient(y), pol.monomial_coefficient(x*y)]

    assert gcd(pol.coefficients()) == 1
    W = abs(coeffs[-1])*X*Y
    assert gcd(coeffs[0],X*Y) == 1
    u = W + ((1-W) % abs(coeffs[0]))
    N = u*(X*Y)^k
    assert gcd(N,coeffs[0]) == 1
    print((2*int(N^(3/4))*X).nbits())
    #print(N)
    mult = inverse_mod(coeffs[0],N)
    coeffs2 = [(coeffs[0]*mult) % N, (coeffs[1]*mult) % N, (coeffs[2]*mult) % N, (coeffs[3]*mult) % N]
    assert coeffs2[0] == 1
    polq = coeffs2[0] + coeffs2[1]*x + coeffs2[2]*y + coeffs2[3]*x*y
    polqq = {}
    for i in range(k+1):
        for j in range(k+1):
            polqq[(i,j)] = polq*x^i*y^j*X^(k-i)*Y^(k-j)
    for i in range(k+2):
        for j in range(k+2):
            if 0 <= i <= k and 0 <= j <= k:
                continue
            polqq[(i,j)] = x^i*y^j*N

    monomials = []
    for i in polqq:
        for j in polqq[i].monomials():
            if j not in monomials:
                monomials.append(j)
    monomials.sort()
    #print(monomials)

    L = matrix(ZZ,len(monomials))

    polqqq = list(polqq.values())
    polqqq.sort()
    for i in range(len(monomials)):
        for j in range(len(monomials)):
            L[i,j] = polqqq[i](X*x,Y*y).monomial_coefficient(monomials[j])
    L = matrix(ZZ,sorted(L,reverse=True))

    print(L.apply_map(lambda x: x.nbits()).str())

    print(L.det().nbits())
    L = L.LLL()
    print(L.apply_map(lambda x: x.nbits()).str())

    for ind in range(L.nrows()):
        PR.<w,z> = PolynomialRing(ZZ)
        pol1 = pol(w,z)
        pol2 = PR(sum((i*j) for i,j in zip(L[ind],monomials))(w/X,z/Y))
        #print(pol2,pol1)
        #print(pol1.parent(), pol2.parent())
        PR.<q> = PolynomialRing(ZZ)
        rr = pol1.resultant(pol2)
        #print(rr)
        if rr.is_zero() or rr.monomials() == [1]:
            print "not independent"
            continue
        rr = rr(q, q)
        #print(rr)
        if len(rr.roots()) > 0:
            print "!!!", rr.roots()
            print(pol(q,rr.roots()[0][0]).roots()[0][0])
            print(rr.roots()[0][0])

def main():
    # Example: recover p,q prime given n=pq and the least significant bits of p

    nbits = 512 # bitlength of primes
    p = random_prime(2^nbits-1, proof=False, lbound=2^(nbits-1))
    q = random_prime(2^nbits-1, proof=False, lbound=2^(nbits-1))
    n = p*q

    lbits = 300 # lower bits of p
    ln = 2^lbits
    p0 = p % ln

    x0 = p // ln # upper bits of p
    y0 = q // ln # upper bits of q
    print 'x0:', x0
    print 'y0:', y0

    # Recovery starts here
    q0 = (n * inverse_mod(p0,ln)) % ln
    assert q0 == q % ln
    X = Y = next_prime(2^(nbits+1-lbits), proof=False) # bounds on x0 and y0

    P.<x,y> = PolynomialRing(ZZ)
    pol = (ln*x+p0)*(ln*y+q0) - n # Should have a root at (x0,y0)

    coron(pol, X, Y, k=2)

if __name__ == '__main__':
    main()
