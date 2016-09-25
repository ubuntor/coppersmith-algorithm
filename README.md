# Coppersmith's Algorithm
Implements Coron's reformulation of Coppersmith's algorithm for finding small integer roots of a bivariate polynomial modulo an integer.

Paper: http://www.jscoron.fr/publications/bivariate.pdf

Used in CSAW CTF Quals 2016 to solve Still Broken Box. (BTW, if you want an implementation of a crypto algorithm, write a crypto CTF challenge that needs it and read writeups.)

**Warning: CTF Quality Code!**

## Why?
Why not?
## Doesn't Sage provide this with `small_roots()`?
`small_roots()` only works with univariate polynomials. (which still would have saved me a lot of time in CSAW...)
