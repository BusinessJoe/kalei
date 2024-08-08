use num_rational::Rational64;

use crate::{matrix::Matrix, poly::{Polynomial, RootLocation}};

#[derive(Debug, Clone)]
pub struct Algebraic {
    poly: Polynomial,
    location: RootLocation,
}

impl std::ops::Mul for Algebraic {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mat1 = Matrix::companion(self.poly);
        let mat2 = Matrix::companion(rhs.poly);

        let mat_prod = mat1.kronecker_product(mat2);

        let poly = mat_prod.char_poly();

        todo!();
    }
}

pub fn sqrt(x: Rational64) -> Algebraic {
    let poly = Polynomial::new(&[-x, (0).into(), (1).into()]);
    let location = RootLocation::Interval((0).into(), x);

    Algebraic {
        poly,
        location,
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multiply() {
        let root2 = sqrt((2).into());
        let _prod = root2.clone() * root2;
    }
}
