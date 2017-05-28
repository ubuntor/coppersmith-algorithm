## Common

This is an instance of common prime RSA, which has a large common factor
between `p-1` and `q-1` to allow for a small private exponent without being
vulnerable to Wiener's attack.

The intended solution is to implement Jochemsz and May's attack in
"A Strategy for Finding Roots of Multivariate Polynomials with New Applications
in attacking RSA Variants".

Players first have to find the paper of the attack - there are three potential
papers when googling, and two of them can (hopefully) be filtered out by the
bounds given in `generate.sage`. If needed, a hint can be given out to filter
out those papers. The correct paper is freely available as a PDF, although
some google-fu may be required.

Players then have to implement the trivariate version of Coppersmith's
algorithm. The reference solution implements an extended version of Coron's
simplification of the bivariate version of Coppersmith's algorithm.

This challenge will be **tough**. It's around the same level of difficulty as
BKP2017's SIDH chall, albeit a bit easier due to the implementation being
easier.
