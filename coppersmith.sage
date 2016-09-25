n = 123541066875660402939610015253549618669091153006444623444081648798612931426804474097249983622908131771026653322601466480170685973651622700515979315988600405563682920330486664845273165214922371767569956347920192959023447480720231820595590003596802409832935911909527048717061219934819426128006895966231433690709
r = nn = 2**300

pp = 8977715745614186723926722135445355667945980810121704263876368391084600429576492555888208633462769980430744297246292311377570428122059389430802674588816921
qq = 10505766383528426528268767841990920439414390641257609651194038542700275920176017420270276503125184671266186854248085247310949691722212477821478080650137921
n = pp*qq
p0 = pp % nn
x0 = pp // nn
y0 = qq // nn
ps = [p0]

for p0 in ps:
    q0 = (n * inverse_mod(p0,nn)) % nn
    print(p0,q0)
    #assert q0 == qq % nn
    X = Y = next_prime(2**(512+1-300), proof=False)
    P.<x,y> = PolynomialRing(ZZ)
    # 0 = (r*x+p0)*(r*y+q0) - n
    # 0 = (p0*q0-n) + (q0*r)*x + (p0*r)*y + (r*r)*x*y
    coeffs = [(p0*q0-n)//r, (q0*r)//r, (p0*r)//r, (r*r)//r]
    polp = coeffs[0] + coeffs[1]*x + coeffs[2]*y + coeffs[3]*x*y
    assert gcd(coeffs) == 1
    W = abs(coeffs[-1])*X*Y
    assert gcd(coeffs[0],X*Y) == 1
    k = 2 # increase if it fails
    u = W + ((1-W) % abs(coeffs[0]))
    N = u*(X*Y)^k
    assert gcd(N,coeffs[0]) == 1
    print(len(str(2*int(N^(3/4))*X)))
    #print(N)
    mult = inverse_mod(coeffs[0],N)
    coeffs2 = [(coeffs[0]*mult) % N, (coeffs[1]*mult) % N, (coeffs[2]*mult) % N, (coeffs[3]*mult) % N]
    assert coeffs2[0] == 1
    polq = coeffs2[0] + coeffs2[1]*x + coeffs2[2]*y + coeffs2[3]*x*y
    #assert polq(x0,y0) % N == 0
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

    print(matrix([[len(str(j)) for j in i] for i in L]))

    '''
    for i in range(L.nrows()):
        poltest = sum((i*j) for i,j in zip(L[i],monomials))
        print(poltest)
        assert poltest(x0/X,y0/Y) % N == 0
    print('matrix good?')
    '''

    print(len(str(L.det())))
    L = L.LLL()
    print(matrix([[len(str(j)) for j in i] for i in L]))

    for ind in range(L.nrows()):
        PR.<w,z> = PolynomialRing(ZZ)
        pol1 = coeffs[0] + coeffs[1]*w + coeffs[2]*z + coeffs[3]*w*z
        pol2 = PR(sum((i*j) for i,j in zip(L[ind],monomials))(w/X,z/Y))
        #print(pol2,pol1)
        #print(pol1.parent(), pol2.parent())
        PR.<q> = PolynomialRing(ZZ)
        #assert pol1(x0,y0) % N == 0
        #assert pol2(x0,y0) % N == 0
        rr = pol1.resultant(pol2)
        #print(rr)
        #if rr(x0,y0) != 0:
        #    print "not good"
        #    continue
        if rr.is_zero() or rr.monomials() == [1]:
            print "not independent"
            continue
        rr = rr(q, q)
        #print(rr)
        if len(rr.roots()) > 0:
            print "!!!", rr.roots()
    print
