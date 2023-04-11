def coron(pol, X, Y, M=None, k=2, debug=False):
    """
    Returns all small roots of pol over the integers (modulo M if it is given).

    Applies Coron's reformulation of Coppersmith's algorithm for finding small
    integer roots of bivariate polynomials modulo an integer.

    Args:
        pol: The polynomial to find small integer roots of.
        X: Upper limit on x.
        Y: Upper limit on y.
        M: Modulus. If M==None, then pol is considered over the integers.
        k: Determines size of lattice. Increase if the algorithm fails.
        debug: Turn on for debug print stuff.

    Returns:
        A list of successfully found roots [(x0,y0), ...].

    Raises:
        ValueError: If pol is not bivariate
    """

    if pol.nvariables() != 2:
        raise ValueError("pol is not bivariate")

    P.<x,y> = PolynomialRing(ZZ)
    pol = P( pol(x,y) )

    # removing common factor of the coefficients
    if M:
        M //= gcd(M,pol.content())
    pol //= pol.content()

    if len(pol.factor()) > 1:
        raise ValueError("pol is reducible")

    # Handle case where pol(0,0) == 0
    xoffset = 0
    while pol(xoffset,0) == 0 or (M and gcd(pol(xoffset,0),M) != 1):
        xoffset += 1
    if debug:
        print("Offset:", xoffset)

    pol = pol(x+xoffset,y)
    p00 = pol(0,0)

    # Handle case where gcd(pol(0,0),X*Y) != 1
    while gcd(p00, X) != 1:
        X = next_prime(X, proof=False)

    while gcd(p00, Y) != 1:
        Y = next_prime(Y, proof=False)

    delta = max(pol.degree(x),pol.degree(y))            # maximum degree of any variable

    if M:
        u = M
    else:
        W = max(abs(i) for i in pol(x*X,y*Y).coefficients())
        u = W + ((1-W) % abs(p00))

    N = u*(X*Y)^k                                       # modulus for polynomials

    # Construct polynomials
    p00inv = inverse_mod(p00,N)
    polq = P( sum((i*p00inv % N)*j for i,j in pol) )

    polynomials = [ polq * x^i * y^j * X^(k-i) * Y^(k-j) if (0 <= i <= k and 0 <= j <= k) 
                    else x^i * y^j * N for i in range(delta+k+1) for j in range(delta+k+1) ]

    # Make list of monomials for matrix indices
    monomials = sorted( set( sum( (i.monomials() for i in polynomials), [] ) ) )

    # Construct lattice spanned by polynomials with xX and yY
    L = matrix(ZZ, len(monomials), lambda i,j: polynomials[i](X*x,Y*y).monomial_coefficient(monomials[j]) )

    # makes lattice upper triangular
    # probably not needed, but it makes debug output pretty
    L = matrix(ZZ,sorted(L,reverse=True))

    if debug:
        print("Bitlengths of matrix elements (before reduction):")
        print(L.apply_map(lambda x: x.nbits()).str())

    L = L.LLL()

    if debug:
        print("Bitlengths of matrix elements (after reduction):")
        print(L.apply_map(lambda x: x.nbits()).str())

    roots = []

    for i in range(L.nrows()):
        if debug:
            print(f"Trying row {i}")

        # i'th row converted to polynomial dividing out X and Y
        pol2 = P(sum(map(mul, zip(L[i],monomials)))(x/X,y/Y))

        r = pol.resultant(pol2, y)

        if r.is_constant():                             # not independent
            continue

        for x0, _ in r.univariate_polynomial().roots():
            if x0+xoffset in [i[0] for i in roots]:
                continue
            if debug:
                print( "Potential x0:",x0 )
            for y0, _ in pol(x0,y).univariate_polynomial().roots():
                if debug:
                    print( "Potential y0:",y0 )
                v = pol(x0,y0)
                if M:
                    v %= M
                if v == 0 and (x0+xoffset,y0) not in roots:
                    roots.append((x0+xoffset,y0))
    return roots

def main():
    # Example 1: recover p,q prime given n=pq and the lower bits of p
    print("---EXAMPLE 1---")

    nbits = 512 # bitlength of primes
    p = random_prime(2^nbits-1, proof=False, lbound=2^(nbits-1))
    q = random_prime(2^nbits-1, proof=False, lbound=2^(nbits-1))
    n = p*q

    lbits = 300 # number of lower bits of p
    ln = 2^lbits
    p0 = p % ln

    x0 = p // ln # upper bits of p
    y0 = q // ln # upper bits of q

    print('p =',p)
    print('q =',q)
    print('x0 =',x0)
    print('y0 =',y0)

    print()
    print('Given:')
    print('n =',n)
    print('p0 =',p0)

    # Recovery starts here
    q0 = (n * inverse_mod(p0,ln)) % ln
    assert q0 == q % ln
    X = Y = 2^(nbits+1-lbits) # bounds on x0 and y0

    P.<x,y> = PolynomialRing(ZZ)
    pol = (ln*x+p0)*(ln*y+q0) - n # Should have a root at (x0,y0)

    x0_2, y0_2 = coron(pol, X, Y, k=2, debug=True)[0]
    p_2 = x0_2*ln + p0
    q_2 = y0_2*ln + q0

    print()
    print('Recovered:')
    print('x0 =',x0_2)
    print('y0 =',y0_2)
    print('p =',p_2)
    print('q =',q_2)

    # Example 2: recover p,q prime given n=pq and the upper bits of p
    # This can be done with a univariate polynomial and Howgrave-Graham,
    # but this is another way to do it with a bivariate polynomial.
    print("---EXAMPLE 2---")

    nbits = 512 # bitlength of primes
    p = random_prime(2^nbits-1, proof=False, lbound=2^(nbits-1))
    q = random_prime(2^nbits-1, proof=False, lbound=2^(nbits-1))
    n = p*q

    lbits = (512-300) # number of masked bits of p
    ln = 2^lbits
    p0 = p // ln

    x0 = p % ln # lower bits of p
    y0 = q % ln # lower bits of q

    print('p =',p)
    print('q =',q)

    print()
    print('Given:')
    print('n =',n)
    print('p0 =',p0)

    # Recovery starts here
    q0 = floor(n / (p0*ln))//ln
    X = Y = 2^(lbits+2) # bounds on x0 and y0
    P.<x,y> = PolynomialRing(ZZ)
    # Should have a root at (x0,y0) +/- some bits of q0
    pol = (x+p0*ln)*(y+q0*ln) - n

    x0_2, y0_2 = coron(pol, X, Y, k=2, debug=True)[0]
    p_2 = p0*ln + x0_2
    q_2 = q0*ln + y0_2

    print()
    print('Recovered:')
    print('x0 =',x0_2)
    print('y0 =',y0_2)
    print('p =',p_2)
    print('q =',q_2)

if __name__ == '__main__':
    main()
