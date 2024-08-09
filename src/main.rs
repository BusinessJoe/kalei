use kalei::ipoly::{IPolynomial, euclid_gcd};

fn main() {
    let a = IPolynomial::new(&[-9, 0, 2]);
    let b = IPolynomial::new(&[0, 4]);
    let gcd = euclid_gcd(&a, &b);

    dbg!(gcd);
}
