mod algebraic;
mod poly;

use poly::Polynomial;

fn main() {
    let p = Polynomial::new(&[(1).into(), (-4).into(), (-4).into(), (8).into()]);

    for loc in p.real_root_locations() {
        println!("{:?}", loc);
    }
}
