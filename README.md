# Montgomery Elliptic Curve Cryptography (mont-ecc)

## About this Repository

This repository offers a comprehensive implementation of elliptic curves, focusing on two primary forms:
- **Short Weierstrass Curve**: Utilizes Jacobian coordinates to enhance computational efficiency.
- **Twisted Edwards Curves**: Employs Extended coordinates, catering to a broader range of curve definitions.

The formulas used for both curve types adhere to the standards set by the Explicit Formulas Database ([EFD](https://hyperelliptic.org/EFD/index.html)).

## Example Usages

The following examples illustrate basic operations such as point addition, doubling, and scaling on different types of curves:

```python
# Example 1: Short Weierstrass Curve
p = 29  # Define the finite field F = Z/pZ
coeffs = [2, 20]  # Curve coefficients

# Initializing the domain with specific optimization settings
domain = Montgomery.factory(mod=p, mul_opt='real0').build()
curve = ShortWeierstrassCurve(domain, coeffs)

# Example 2: Twisted Edwards Curve
coeffs = [13, 21]
edward = TwistedEdwardsCurve(domain, coeffs)

# Perform curve operations such as point addition, doubling, and scalar multiplication

# <skipped> refer to unit tests for details
```

## Notes:

For the Twisted Edwards curves, this implementation defaults to the "unified" addition and doubling formulas. This choice allows the library to support more generalized forms of Twisted Edwards curves, extending its applicability beyond standard protocol implementations. For standard curve operations, the dedicated formulas can also be used without issues as detailed in Reference 4.

## References
1. [NIST SP 800-186](https://csrc.nist.gov/pubs/sp/800/186/final)
2. [RFC7748](https://www.rfc-editor.org/rfc/rfc7748)
3. [Twisted Edwards Curves](https://eprint.iacr.org/2008/013.pdf)
4. [Twisted Edwards Curves Revisited](https://eprint.iacr.org/2008/522.pdf)
