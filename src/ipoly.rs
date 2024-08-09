use num_integer::{gcd, lcm};
use num_rational::Rational64;
use itertools::{Itertools, ZipLongest};

type Integer = i64;

#[derive(Debug, Clone, PartialEq, Eq)]
/// Polynomial over the integers
struct IPolynomial {
    coefs: Vec<Integer>,
}

impl IPolynomial {
    pub fn new(coefs: &[Integer]) -> Self {
        Self {
            coefs: coefs.to_vec(),
        }
    }

    pub fn degree(&self) -> usize {
        // We always have a positive number of coefficients since a polynomial
        // always has at least one coefficient (even if it's just [0]).
        self.coefs.len() - 1
    }

    pub fn leading_coef(&self) -> Integer {
        // We always have a positive number of coefficients since a polynomial
        // always has at least one coefficient (even if it's just [0]).
        *self.coefs.last().unwrap()
    }

    pub fn derivative(&self) -> Self {
        let f_prime_coefs: Vec<i64> = self.coefs[1..]
            .iter()
            .enumerate()
            .map(|(i, c)| i64::try_from(i + 1).unwrap() * c)
            .collect();
        IPolynomial::new(&f_prime_coefs)
    }
}

impl std::ops::Sub for IPolynomial {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut coefs: Vec<_> = self.coefs.iter().zip_longest(rhs.coefs).map(|zip| {
            match zip {
                itertools::EitherOrBoth::Both(lc, rc) => lc - rc,
                itertools::EitherOrBoth::Left(lc) => *lc,
                itertools::EitherOrBoth::Right(rc) => -rc,
            }
        }).collect();

        while coefs.len() > 1 && Some(&0) == coefs.last() {
            coefs.pop();
        }

        Self::new(&coefs)
    }
}

// TODO: write proper citation
/// Implemented from:
/// A Verified Implementation of Algebraic Numbers In Isabelle/HOL
struct Algebraic {
    pub poly: IPolynomial,
}

impl Algebraic {
    pub fn from_rational(mut numer: Integer, mut denom: Integer) -> Option<Self> {
        if denom == 0 {
            return None;
        }

        // Remove common factors
        let divisor = gcd(numer, denom);
        numer /= divisor;
        denom /= divisor;

        // Ensure denominator is positive
        if denom < 0 {
            numer = -numer;
            denom = -denom;
        }

        // Now the number can be represented by dx - n
        let poly = IPolynomial::new(&[-numer, denom]);
        Some(Self { poly })
    }

    pub fn from_integer(x: Integer) -> Self {
        // from_rational only fails if the denominator is zero
        Self::from_rational(x, 1).unwrap()
    }
}

impl std::ops::Neg for Algebraic {
    type Output = Self;

    fn neg(self) -> Self::Output {
        // Make the substitution x -> -x
        let mut poly = self.poly.clone();
        let degree = poly.degree();
        for (idx, coef) in poly.coefs.iter_mut().enumerate() {
            if idx % 2 == 1 {
                *coef *= -1
            }
            if degree % 2 == 1 {
                *coef *= -1
            }
        }
        Self { poly }
    }
}

/// Returns polynomials (q, r) such that a = qb + r.
/// This algorithm only works for polynomials over fields, so we return
/// rational coefficients instead of integral ones.
fn euclid_div_rational(a: &[Rational64], b: &[Rational64]) -> (Vec<Rational64>, Vec<Rational64>) {
    let mut q: Vec<Rational64> = vec![(0).into(); a.len()];
    let mut r: Vec<Rational64> = a.iter().map(|&c| c.into()).collect();
    let d = b.len() - 1;
    let c = b.last().unwrap();

    while r.len() - 1 >= d {
        let s = r.last().unwrap() / c;
        let s_exp = r.len() - 1 - d;
        q[s_exp] += s;

        // r = r - sb
        for (idx, &coef) in b.iter().enumerate() {
            r[idx + s_exp] -= Rational64::from(coef) * s
        }

        // Remove leading zeros from r
        while r.len() > 1 && Some(&Rational64::ZERO) == r.last() {
            r.pop();
        }
    }

    // Remove leading zeros from q
    while q.len() > 1 && Some(&Rational64::ZERO) == q.last() {
        q.pop();
    }

    (q, r)
}

/// Assumes a is divisible by b
fn ipoly_div(a: &IPolynomial, b: &IPolynomial) -> IPolynomial {
    let a_rat: Vec<_> = a.coefs.iter().map(|&c| Rational64::from_integer(c)).collect();
    let b_rat: Vec<_> = b.coefs.iter().map(|&c| Rational64::from_integer(c)).collect();

    // Assume that the remainder is 0
    let (q, _) = euclid_div_rational(&a_rat, &b_rat);

    let coefs: Vec<_> = q.into_iter().map(|c| {
        debug_assert!(c.is_integer());
        *c.numer()
    }).collect();

    IPolynomial::new(&coefs)
}

fn euclid_gcd_rational(a: &[Rational64], b: &[Rational64]) -> Vec<Rational64> {
    if a.len() < b.len() {
        return euclid_gcd_rational(b, a);
    }
    if b == [Rational64::ZERO] {
        return a.to_vec();
    }

    let (_, r) = euclid_div_rational(a, b);
    euclid_gcd_rational(b, &r)
}

/// Returns the GCD of two polynomials.
/// Implemented via Euclidean algorithm.
fn euclid_gcd(a: &IPolynomial, b: &IPolynomial) -> IPolynomial {
    // Extend to rationals for gcd algorithm
    let a_rat: Vec<_> = a.coefs.iter().map(|&c| Rational64::from_integer(c)).collect();
    let b_rat: Vec<_> = b.coefs.iter().map(|&c| Rational64::from_integer(c)).collect();

    let gcd_rat = euclid_gcd_rational(&a_rat, &b_rat);

    dbg!(&gcd_rat);

    // Find lcm of denominators
    let mut denom_lcm = 1;
    for c in &gcd_rat {
        denom_lcm = lcm(denom_lcm, *c.denom());
    }

    let mut coefs: Vec<_> = gcd_rat.into_iter().map(|mut coef| {
        coef *= denom_lcm;
        debug_assert!(coef.is_integer());
        *coef.numer()
    }).collect();

    // Divide coefs by gcd
    let mut coef_gcd = coefs[0];
    for c in &coefs {
        coef_gcd = gcd(coef_gcd, *c);
    }

    coefs.iter_mut().for_each(|c| *c /= coef_gcd);
    IPolynomial::new(&coefs)
}


/// Returns the square free decomposition of the given polynomial.
/// Implemented via Yun's algorithm.
fn square_free_decomposition(f: &IPolynomial) -> Vec<IPolynomial> {
    let f_prime = f.derivative();

    let mut outputs: Vec<IPolynomial> = vec![
        euclid_gcd(&f, &f_prime)
    ];

    let mut b = ipoly_div(&f, &outputs[0]);
    let mut c = ipoly_div(&f_prime, &outputs[0]);
    let mut d = c - b.derivative();
    let mut iter = 0;

    while b.coefs != [1] {
        dbg!(&b, &d);
        panic!("bruh");
        let a = euclid_gcd(&b, &d);
        b = ipoly_div(&b, &a);
        c = ipoly_div(&d, &a);
        d = c - b.derivative();
        outputs.push(a);
        iter += 1;
    }

    outputs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_euclid_div() {
        // Test (x^2 + 4) / (x + 1)
        // Expect x^2 + 4 = (x - 1)(x + 1) + 5
        let a = [4, 0, 1].map(Rational64::from_integer);
        let b = [1, 1].map(Rational64::from_integer);

        let (q, r) = euclid_div_rational(&a, &b);

        dbg!(&q, &r);

        assert_eq!(q, [(-1).into(), (1).into()]);
        assert_eq!(r, [(5).into()]);
    }

    #[test]
    fn test_euclid_div_leading_coef_not_one() {
        // Test (2x^2 + 4) / (2x + 1)
        // Expect 2x^2 + 4 = (x - 1/2)(2x + 1) + 9/2
        let a = [4, 0, 2].map(Rational64::from_integer);
        let b = [1, 2].map(Rational64::from_integer);

        let (q, r) = euclid_div_rational(&a, &b);
        dbg!(&q, &r);

        assert_eq!(q, [Rational64::new(-1, 2), (1).into()]);
        assert_eq!(r, [Rational64::new(9, 2)]);
    }

    #[test]
    fn test_euclid_gcd_simple() {
        // Find gcd(x, 1)
        // Expect 1
        let a = IPolynomial::new(&[0, 1]);
        let b = IPolynomial::new(&[1]);

        let gcd = euclid_gcd(&a, &b);
        dbg!(&gcd);
        
        assert_eq!(gcd.coefs, [1]);
    }

    #[test]
    fn test_euclid_gcd() {
        // Find gcd(x^2 + 7x + 6, x^2 - 5x - 6)
        // Expect x + 1
        let a = IPolynomial::new(&[6, 7, 1]);
        let b = IPolynomial::new(&[-6, -5, 1]);

        let gcd = euclid_gcd(&a, &b);
        dbg!(&gcd);
        
        assert_eq!(gcd.coefs, [1, 1]);
    }

    #[test]
    fn test_euclid_gcd_quartic() {
        // Find gcd(x^4 - 4x^3 + 4x^2 - 3x + 14, x^4 + 8x^3 + 12x^2 + 17x + 6)
        // Expect x^2 + x + 2
        let a = IPolynomial::new(&[14, -3, 4, -4, 1]);
        let b = IPolynomial::new(&[6, 17, 12, 8, 1]);

        let gcd = euclid_gcd(&a, &b);
        dbg!(&gcd);
        
        assert_eq!(gcd.coefs, [2, 1, 1]);
    }

    #[test]
    fn test_derivative() {
        // x^2 -> 2x
        let a = IPolynomial::new(&[0, 0, 1]);

        assert_eq!(a.derivative(), IPolynomial::new(&[0, 2]))
    }

    #[test]
    fn test_square_free_decomposition() {
        // x^2 -> x
        let a = IPolynomial::new(&[0, 0, 1]);
        let square_free = square_free_decomposition(&a);
        assert_eq!(square_free, vec![IPolynomial::new(&[0, 1])])
    }
}
