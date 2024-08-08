use num_integer::gcd;

type Integer = i64;

#[derive(Debug, Clone)]
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
        // TODO: does this work for polynomials of odd degree? won't they have
        // a negative leading coefficient after this?
        let mut poly = self.poly.clone();
        for (idx, coef) in poly.coefs.iter_mut().enumerate() {
            if idx % 2 == 1 {
                *coef *= -1
            }
        }
        Self { poly }
    }
}
