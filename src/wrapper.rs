use super::convex::{ConvexHull, ErrorKind};
use super::util::det;
use num_bigint::{BigInt, ToBigInt};
use num_traits::Float;

#[derive(Clone, Debug)]
// A wrapper that holds a given points as a bigint.
// This can compute a convex hull robustly, although the computation is slower.
pub struct ConvexHullWrapper<T: Float + ToBigInt> {
    pub inner: ConvexHull<BigInt>,
    pub conversion_factor: T,
}

impl<T: Float + ToBigInt> ConvexHullWrapper<T> {
    pub fn try_new(points: &[Vec<T>], max_iter: Option<usize>) -> Result<Self, ErrorKind> {
        let conversion_factor =T::epsilon();
        let int_points:Vec<_> = to_bigint_points(points, conversion_factor);
        let convex_hull = ConvexHull::try_new(&int_points, 0, max_iter)?;
        Ok(
            Self{
                inner: convex_hull,
                conversion_factor,
            }
        )
    }
    pub fn add_points(
        &mut self,
        points: &[Vec<T>],
    ) -> Result<(), ErrorKind>{
        let int_points = to_bigint_points(points, self.conversion_factor);
        self.inner.add_points(&int_points, 0)
    }
    pub fn vertices_indices(&self) -> (Vec<Vec<T>>, Vec<usize>) {
        let (int_v, i) = self.inner.vertices_indices();
        let v = to_float_points(&int_v, self.conversion_factor);
        (v,i)
    }
    pub fn volume(&self) -> T {
        let dim = self.inner.points[0].len();
        let (c_hull_vertices, c_hull_indices) = self.vertices_indices();
        let mut reference_point = c_hull_vertices[c_hull_indices[0]].to_vec();
        reference_point.push(T::one());
        let mut volume = T::zero();
        for i in (dim..c_hull_indices.len()).step_by(dim) {
            let mut mat = Vec::new();
            for j in 0..dim {
                let mut row = c_hull_vertices[c_hull_indices[i + j]].to_vec();
                row.push(T::one());
                mat.push(row);
            }
            mat.push(reference_point.to_vec());
            volume = volume + det(&mat);
        }
        let factorial = {
            let mut result = T::one();
            let mut m = T::one() + T::one();
            let mut n = dim;
            while n > 1 {
                result = result * m.clone();
                n = n - 1;
                m = m + T::one();
            }
            result
        };
        volume / factorial
    }
    pub fn support_point(&self, direction: &[T]) -> Result<Vec<T>, ErrorKind> {
        let d = to_bigint_points(&[direction.to_vec()], self.conversion_factor);
        let p = to_float_points(&[self.inner.support_point(&d[0])?], self.conversion_factor);
        Ok(p[0].to_vec())
    }
}

fn to_bigint_points<T: Float+ToBigInt>(points: &[Vec<T>], conversion_factor: T) ->Vec<Vec<BigInt>>{
    let mut int_points = Vec::new();
    for point in points {
        int_points.push(
            point
                .iter()
                .map(|&p| p / conversion_factor)
                .map(|p| p.to_bigint().expect("cannot convert to bigint"))
                .collect::<Vec<_>>(),
        );
    }
    int_points
}

fn to_float_points<T: Float+ToBigInt>(points: &[Vec<BigInt>], conversion_factor: T) ->Vec<Vec<T>>{
    let mut float_points = Vec::new();
    for point in points {
        float_points.push(
            point
                .iter()
                .map(|p| T::from(p.clone()).expect("cannot convert to bigint") * conversion_factor)
                .collect::<Vec<_>>(),
        );
    }
    float_points
}

#[test]
fn wrapper_cube_test() {
    let p1 = vec![1.0, 1.0, 1.0];
    let p2 = vec![1.0, 1.0, -1.0];
    let p3 = vec![1.0, -1.0, 1.0];
    let p4 = vec![1.0, -1.0, -1.0];
    let p5 = vec![-1.0, 1.0, 1.0];
    let p6 = vec![-1.0, 1.0, -1.0];
    let p7 = vec![-1.0, -1.0, 1.0];
    let p8 = vec![-1.0, -1.0, -1.0];
    let (_v, i) = ConvexHullWrapper::try_new(&[p1, p2, p3, p4, p5, p6, p7, p8], None)
        .unwrap()
        .vertices_indices();
    assert_eq!(i.len(), 6 * 2 * 3);
}

#[test]
fn wrapper_cube_volume_test() {
    let p1 = vec![2.0, 2.0, 2.0];
    let p2 = vec![2.0, 2.0, 0.0];
    let p3 = vec![2.0, 0.0, 2.0];
    let p4 = vec![2.0, 0.0, 0.0];
    let p5 = vec![0.0, 2.0, 2.0];
    let p6 = vec![0.0, 2.0, 0.0];
    let p7 = vec![0.0, 0.0, 2.0];
    let p8 = vec![0.0, 0.0, 0.0];
    let cube = ConvexHullWrapper::try_new(&[p1, p2, p3, p4, p5, p6, p7, p8], None).unwrap();
    assert_eq!(cube.volume(), 8.0);
}
