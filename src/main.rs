mod algebraic;
mod poly;

use num_rational::Rational64;
use poly::Polynomial;

fn main() {
    // let p = Polynomial::new(&[
    //     Rational64::from_integer(-2),
    //     Rational64::ZERO,
    //     Rational64::from_integer(1),
    // ]);
    let p = Polynomial::new(&[
        Rational64::from_integer(-1),
        Rational64::from_integer(-4),
        Rational64::from_integer(4),
        Rational64::from_integer(8),
    ]);

    for loc in p.real_root_locations() {
        println!("{:?}", loc);
    }
}
