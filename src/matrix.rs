use num_rational::Rational64;

use crate::poly::Polynomial;

#[derive(Debug)]
pub struct Matrix<T> {
    rows: usize,
    columns: usize,
    elements: Vec<T>,
}

impl<T> Matrix<T> {
    pub fn new(rows: usize, columns: usize, elements: Vec<T>) -> Self {
        if rows * columns != elements.len() {
            panic!("wrong number of elements")
        }

        Self {
            rows,
            columns,
            elements,
        }
    }

    pub fn get(&self, row: usize, column: usize) -> &T {
        let index = row * self.columns + column;
        &self.elements[index]
    }

    pub fn get_mut(&mut self, row: usize, column: usize) -> &mut T {
        let index = row * self.columns + column;
        &mut self.elements[index]
    }

    pub fn mat_mul(&self, right: &Self) -> Self {
        todo!()
    }
}

impl Matrix<Rational64> {
    /// Return the companion matrix of the given polynomial
    pub fn companion(poly: Polynomial) -> Self {
        // Make polynomial monic
        let mut monic = poly.clone();
        let leading_coef = monic.coefs[monic.coefs.len() - 1];
        for coef in monic.coefs.iter_mut() {
            *coef /= leading_coef;
        }

        // The companion matrix of P is a square matrix with side length equal to the degree of P.
        let rows = monic.degree();
        let elements: Vec<Rational64> = vec![(0).into(); rows * rows];

        // Assign 1s along diagonal
        let mut matrix = Matrix::new(rows, rows, elements);
        for c in 0..rows - 1 {
            let r = c + 1;
            *matrix.get_mut(r, c) = (1).into();
        }

        // Assign coefs along rightmost column
        for r in 0..rows {
            let coef = monic.coefs[r];
            *matrix.get_mut(r, rows - 1) = coef;
        }

        matrix
    }

    pub fn kronecker_product(&self, other: Self) -> Self {
        let rows = self.rows * other.rows;
        let columns = self.columns * other.columns;
        let elements: Vec<Rational64> = vec![(0).into(); rows * columns];
        let mut matrix = Matrix::new(rows, columns, elements);

        for a_row in 0..self.rows {
            let row_start = a_row * other.rows;
            for a_col in 0..self.columns {
                let col_start = a_col * other.columns;

                for b_row in 0..other.rows {
                    let row = row_start + b_row;

                    for b_col in 0..other.columns {
                        let col = col_start + b_col;

                        let a_element = self.get(a_row, a_col);
                        let b_element = other.get(b_row, b_col);

                        *matrix.get_mut(row, col) = a_element * b_element;
                    }
                }
            }
        }

        matrix
    }

    pub fn char_poly(&self) -> Polynomial {
        println!("{:?}", self);
        assert_eq!(self.rows, self.columns, "matrix must be square");
        todo!()
    }
}
