use std::{cmp::max, fmt::Display};

use num_integer::binomial;
use num_rational::Rational64;
use num_traits::Signed;

#[derive(Debug, Clone)]
pub struct Polynomial {
    pub coefs: Vec<Rational64>,
}

#[derive(Debug, Clone)]
pub enum RootLocation {
    Exact(Rational64),
    Interval(Rational64, Rational64),
}

impl RootLocation {
    pub fn contains(&self, value: Rational64) -> bool {
        match self {
            &RootLocation::Exact(x) => value == x,
            &RootLocation::Interval(l, u) => l < value && value < u,
        }
    }
}

impl Display for RootLocation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RootLocation::Exact(x) => write!(f, "{}", x),
            RootLocation::Interval(l, u) => write!(f, "({}, {})", l, u),
        }
    }
}

impl Polynomial {
    pub fn new(coefs: &[Rational64]) -> Self {
        Self {
            coefs: coefs.to_vec(),
        }
    }

    pub fn degree(&self) -> usize {
        self.coefs.len() - 1
    }

    fn count_sign_changes(&self) -> usize {
        let mut count = 0;
        let mut last_sign_pos = None;
        for coef in self.coefs.iter().filter(|&c| c != &Rational64::ZERO) {
            if let Some(last_sign) = last_sign_pos {
                if last_sign != (coef > &Rational64::ZERO) {
                    count += 1;
                }
            }
            last_sign_pos = Some(coef > &Rational64::ZERO);
        }
        count
    }

    /// Makes the substitution x -> x + 1
    fn sub_shift(&self) -> Self {
        let mut coefs = vec![];
        for i in 0..self.coefs.len() {
            let mut coef = Rational64::ZERO;
            for j in i..self.coefs.len() {
                coef += Rational64::from_integer(binomial(j, i).try_into().unwrap()) * self.coefs[j]
            }

            coefs.push(coef)
        }

        Self { coefs }
    }

    /// Makes the substitution x -> 1 / (x + 1)
    fn sub_reciprocal(&self) -> Self {
        let mut coefs = vec![];
        for i in 0..self.coefs.len() {
            let mut coef = Rational64::ZERO;
            for j in i..self.coefs.len() {
                coef += Rational64::from_integer(binomial(j, i).try_into().unwrap())
                    * self.coefs[self.coefs.len() - 1 - j]
            }

            coefs.push(coef)
        }

        Self { coefs }
    }

    /// Implemented with Uspensky's algorithm
    fn find_real_roots_recursive(&self, debug: &str) -> Vec<RootLocation> {
        let sign_changes = self.count_sign_changes();
        if sign_changes == 0 {
            return vec![];
        }
        if sign_changes == 1 {
            let upper = self.root_upper_bound();
            return vec![RootLocation::Interval(0.into(), upper)];
        }

        let mut locations = vec![];

        let a = self.sub_shift();
        if a.coefs[0] == Rational64::ZERO {
            locations.push(RootLocation::Exact(1.into()))
        }

        for transformed_location in a.find_real_roots_recursive(&format!("{debug} x+1")) {
            let location = match transformed_location {
                RootLocation::Exact(x) => RootLocation::Exact(x + 1),
                RootLocation::Interval(l, u) => RootLocation::Interval(l + 1, u + 1),
            };
            locations.push(location);
        }

        let b = self.sub_reciprocal();
        for transformed_location in b.find_real_roots_recursive(&format!("{debug} 1/(x+1)")) {
            let location = match transformed_location {
                RootLocation::Exact(x) => RootLocation::Exact((x + 1).recip()),
                RootLocation::Interval(l, u) => {
                    RootLocation::Interval((u + 1).recip(), (l + 1).recip())
                }
            };
            locations.push(location);
        }

        locations
    }

    /// Factors out roots at x = 0
    fn remove_trivial_zeros(&mut self) -> bool {
        let mut had_zeros = false;
        while self.coefs[0] == Rational64::ZERO {
            self.coefs.remove(0);
            had_zeros = true;
        }
        had_zeros
    }

    /// Returns a rational with greater magnitude than any real zero of this polynomial
    fn root_upper_bound(&self) -> Rational64 {
        let mut bound = Rational64::ZERO;
        for coef in &self.coefs[0..self.coefs.len() - 1] {
            bound = max(coef.abs(), bound);
        }
        bound /= self.coefs.last().unwrap().abs();
        bound + 1
    }

    /// Returns locations which each contain a single real root.
    /// Locations can be an exact value or a lower and upper bound.
    pub fn real_root_locations(&self) -> Vec<RootLocation> {
        let mut locations = vec![];

        // Clone so we can be destructive
        let mut p = self.clone();

        // Factor out trivial root at x = 0
        let had_zero_root = p.remove_trivial_zeros();
        if had_zero_root {
            locations.push(RootLocation::Exact(Rational64::ZERO));
        }

        locations.extend_from_slice(&p.find_real_roots_recursive("pos"));

        // Substitute x -> -x
        for (i, coef) in p.coefs.iter_mut().enumerate() {
            if i % 2 == 1 {
                *coef *= -1;
            }
        }

        locations.extend(
            p.find_real_roots_recursive("neg")
                .into_iter()
                .map(|loc| match loc {
                    RootLocation::Exact(x) => RootLocation::Exact(-x),
                    RootLocation::Interval(l, u) => RootLocation::Interval(-u, -l),
                }),
        );

        locations
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn count_sign_changes() {
        let p = Polynomial::new(&[
            Rational64::from_integer(-2),
            Rational64::ZERO,
            Rational64::from_integer(1),
        ]);
        assert_eq!(p.count_sign_changes(), 1);
    }

    #[test]
    fn count_sign_changes_cos_2pi_over_7() {
        let p = Polynomial::new(&[
            Rational64::from_integer(-1),
            Rational64::from_integer(-4),
            Rational64::from_integer(4),
            Rational64::from_integer(8),
        ]);
        assert_eq!(p.count_sign_changes(), 1);
    }

    #[test]
    fn find_root_of_cos_pi_7_min_poly() {
        let p = Polynomial::new(&[(1).into(), (-4).into(), (-4).into(), (8).into()]);

        let root_locations = p.real_root_locations();
        assert_eq!(root_locations.len(), 3);

        let mut expected_roots = vec![-0.623, 0.223, 0.901];
        for loc in &root_locations {
            let mut root_idx = None;
            for (i, root) in expected_roots.iter().enumerate() {
                if loc.contains(Rational64::approximate_float(*root).unwrap()) {
                    root_idx = Some(i);
                    break;
                }
            }

            if let Some(root_idx) = root_idx {
                expected_roots.remove(root_idx);
            } else {
                panic!("No root found for location {:?}", loc);
            }
        }
    }

    #[test]
    fn root_upper_bound() {
        let p = Polynomial::new(&[(-8).into(), (-16).into(), (-6).into(), (1).into()]);

        // Roots are about -0.7, -1.4, and 8.1
        // Therefore the upper bound should be larger than 8.1
        let bound = p.root_upper_bound();
        assert!(
            bound > 8.into(),
            "bound {} should be larger than {}",
            bound,
            8
        );
    }
}
