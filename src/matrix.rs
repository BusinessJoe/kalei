pub struct Matrix<T> {
    rows: usize,
    columns: usize,
    elements: Vec<T>,
}

impl<T> Matrix<T> {
    pub fn new(rows: usize, columns: usize, elements: Vec<T>) -> Self {
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
