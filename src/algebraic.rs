use num_rational::Rational64;

use crate::poly::Polynomial;


pub struct Algebraic {
    poly: Polynomial,
    lower: Rational64,
    upper: Rational64,
}
