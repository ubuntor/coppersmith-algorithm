def coron(pol, X, Y, k=2, debug=False):
    """
    Returns a small root of pol.

    Applies Coron's reformulation of Coppersmith's algorithm for finding small
    integer roots of bivariate polynomials modulo an integer.

    Args:
        pol: The polynomial to find small integer roots of.
        X: Upper limit on x.
        Y: Upper limit on y.
        k: Determines size of lattice. Increase if the algorithm fails.
        debug: Turn on for debug print stuff.

    Returns:
        (x0,y0) if the algorithm successfully finds a small root (x0,y0).
        (-1,-1) if the algorithm fails.

    Raises:
        ValueError: If pol is not bivariate, pol(0,0) == 0, or
            gcd(pol(0,0),X*Y) != 1.
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

    delta = max(pol.degree(x),pol.degree(y)) # maximum degree of any variable

    p00 = pol(0,0)

    W = max(abs(i) for i in pol(x*X,y*Y).coefficients())
    u = W + ((1-W) % abs(p00))
    N = u*(X*Y)^k # modulus for polynomials

    # Construct polynomials
    p00inv = inverse_mod(p00,N)
    polq = P(sum((i*p00inv % N)*j for i,j in zip(pol.coefficients(),
                                                 pol.monomials())))
    polqq = {}
    for i in range(k+1):
        for j in range(k+1):
            polqq[(i,j)] = polq*x^i*y^j*X^(k-i)*Y^(k-j)
    for i in range(delta+k+1):
        for j in range(delta+k+1):
            if 0 <= i <= k and 0 <= j <= k:
                continue
            polqq[(i,j)] = x^i*y^j*N

    monomials = []
    for i in polqq:
        for j in polqq[i].monomials():
            if j not in monomials:
                monomials.append(j)
    monomials.sort()

    # Construct lattice spanned by polynomials with xX and yY
    L = matrix(ZZ,len(monomials))
    polqqq = list(polqq.values())
    polqqq.sort()
    for i in range(len(monomials)):
        for j in range(len(monomials)):
            L[i,j] = polqqq[i](X*x,Y*y).monomial_coefficient(monomials[j])

    L = matrix(ZZ,sorted(L,reverse=True)) # makes matrix upper triangular

    if debug:
        print "Bitlengths of matrix elements (before reduction):"
        print L.apply_map(lambda x: x.nbits()).str()

    L = L.LLL()

    if debug:
        print "Bitlengths of matrix elements (after reduction):"
        print L.apply_map(lambda x: x.nbits()).str()

    for i in range(L.nrows()):
        if debug:
            print "Trying row %d" % i

        pol1 = pol
        # i'th row converted to polynomial dividing out X and Y
        pol2 = P(sum(map(mul, zip(L[i],monomials)))(x/X,y/Y))

        P2.<z> = PolynomialRing(ZZ)
        r = pol1.resultant(pol2, y)

        if r.is_constant(): # not independent
            continue

        r = r(z,0) # convert to univariate polynomial

        if len(r.roots()) > 0:
            for x0, _ in r.roots():
                if debug:
                    print "Potential x0:",x0
                for y0, _ in P2(pol(x0,z)).roots():
                    if debug:
                        print "Potential y0:",y0
                    if pol(x0,y0) == 0:
                        return (x0,y0)
    return (-1,-1)

def main():
    # Example: recover p,q prime given n=pq and the lower bits of p

    nbits = 512 # bitlength of primes
    p = random_prime(2^nbits-1, proof=False, lbound=2^(nbits-1))
    q = random_prime(2^nbits-1, proof=False, lbound=2^(nbits-1))
    n = p*q

    lbits = 300 # number of lower bits of p
    ln = 2^lbits
    p0 = p % ln

    x0 = p // ln # upper bits of p
    y0 = q // ln # upper bits of q

    print 'p =',p
    print 'q =',q
    print 'x0 =',x0
    print 'y0 =',y0

    print
    print 'Given:'
    print 'n =',n
    print 'p0 =',p0

    # Recovery starts here
    q0 = (n * inverse_mod(p0,ln)) % ln
    assert q0 == q % ln
    X = Y = next_prime(2^(nbits+1-lbits), proof=False) # bounds on x0 and y0

    P.<x,y> = PolynomialRing(ZZ)
    pol = (ln*x+p0)*(ln*y+q0) - n # Should have a root at (x0,y0)

    x0_2, y0_2 = coron(pol, X, Y, k=2, debug=True)
    p_2 = x0_2*ln + p0
    q_2 = y0_2*ln + q0

    print
    print 'Recovered:'
    print 'x0 =',x0_2
    print 'y0 =',y0_2
    print 'p =',p_2
    print 'q =',q_2

if __name__ == '__main__':
    main()
