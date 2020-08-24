use kdtree::distance::squared_euclidean;
use kdtree::KdTree;
use num_traits::{Float, NumOps, Zero};
use std::collections::BTreeMap;
use std::error::Error;

pub(crate) fn min_max_index_each_axis<T: Clone + Zero + PartialOrd>(
    points: &[Vec<T>],
) -> Vec<(usize, usize)> {
    let dim = points[0].len();
    let mut min_index = vec![0; dim];
    let mut max_index = vec![0; dim];
    let mut min = vec![T::zero(); dim];
    let mut max = vec![T::zero(); dim];
    for (index, point) in points.iter().enumerate() {
        for (j, element) in point.iter().enumerate() {
            if index == 0 || *element < min[j] {
                min[j] = element.clone();
                min_index[j] = index;
            }
            if index == 0 || *element > max[j] {
                max[j] = element.clone();
                max_index[j] = index;
            }
        }
    }
    min_index.into_iter().zip(max_index.into_iter()).collect()
}

pub(crate) fn is_same_dimension<T>(points: &[Vec<T>]) -> bool {
    if points.len() == 0 {
        return true;
    }
    let dim = points[0].len();
    if points.iter().skip(1).find(|p| p.len() != dim).is_some() {
        false
    } else {
        true
    }
}

pub(crate) fn det<T: NumOps + Clone>(m: &[Vec<T>]) -> T {
    let row_dim = m.len();
    let column_dim = m[0].len();
    assert_eq!(row_dim, column_dim);
    match column_dim {
        1 => m[0][0].clone(),
        2 => det_2x2(m),
        3 => det_3x3(m),
        4 => det_4x4(m),
        _ => {
            unimplemented!("matrix size should be less than 4 dim to calcurate determinant for now")
        }
    }
}

fn det_2x2<T: NumOps + Clone>(m: &[Vec<T>]) -> T {
    m[0][0].clone() * m[1][1].clone() - m[0][1].clone() * m[1][0].clone()
}

#[rustfmt::skip]
fn det_3x3<T: NumOps+Clone>(m: &[Vec<T>]) -> T {
    m[0][0].clone() * (m[1][1].clone() * m[2][2].clone() - m[1][2].clone() * m[2][1].clone())
  - m[1][0].clone() * (m[0][1].clone() * m[2][2].clone() - m[0][2].clone() * m[2][1].clone())
  + m[2][0].clone() * (m[0][1].clone() * m[1][2].clone() - m[0][2].clone() * m[1][1].clone())
}

#[rustfmt::skip]
fn det_4x4<T: NumOps + Clone>(m: &[Vec<T>]) -> T {
    m[0][0].clone()
        * (   m[1][1].clone() * m[2][2].clone() * m[3][3].clone()
            + m[1][2].clone() * m[2][3].clone() * m[3][1].clone()
            + m[1][3].clone() * m[2][1].clone() * m[3][2].clone()
            - m[1][3].clone() * m[2][2].clone() * m[3][1].clone()
            - m[1][2].clone() * m[2][1].clone() * m[3][3].clone()
            - m[1][1].clone() * m[2][3].clone() * m[3][2].clone())
  - m[1][0].clone()
        * (   m[0][1].clone() * m[2][2].clone() * m[3][3].clone()
            + m[0][2].clone() * m[2][3].clone() * m[3][1].clone()
            + m[0][3].clone() * m[2][1].clone() * m[3][2].clone()
            - m[0][3].clone() * m[2][2].clone() * m[3][1].clone()
            - m[0][2].clone() * m[2][1].clone() * m[3][3].clone()
            - m[0][1].clone() * m[2][3].clone() * m[3][2].clone())
  + m[2][0].clone()
        * (   m[0][1].clone() * m[1][2].clone() * m[3][3].clone()
            + m[0][2].clone() * m[1][3].clone() * m[3][1].clone()
            + m[0][3].clone() * m[1][1].clone() * m[3][2].clone()
            - m[0][3].clone() * m[1][2].clone() * m[3][1].clone()
            - m[0][2].clone() * m[1][1].clone() * m[3][3].clone()
            - m[0][1].clone() * m[1][3].clone() * m[3][2].clone())
  - m[3][0].clone()
        * (   m[0][1].clone() * m[1][2].clone() * m[2][3].clone()
            + m[0][2].clone() * m[1][3].clone() * m[2][1].clone()
            + m[0][3].clone() * m[1][1].clone() * m[2][2].clone()
            - m[0][3].clone() * m[1][2].clone() * m[2][1].clone()
            - m[0][2].clone() * m[1][1].clone() * m[2][3].clone()
            - m[0][1].clone() * m[1][3].clone() * m[2][2].clone())
}

pub fn remove_nearby_points<T>(
    points: &[Vec<T>],
    squared_distance: T,
) -> Result<Vec<Vec<T>>, Box<dyn Error>>
where
    T: Float,
{
    let mut points_map = BTreeMap::new();
    let dim = points[0].len();
    let mut kdtree = KdTree::new(dim);
    for (i, point) in points.iter().enumerate() {
        points_map.insert(i, point);
        kdtree.add(point, i)?;
    }
    let mut live_indices = Vec::new();
    while let Some((id, point)) = points_map.iter().next() {
        let mut remove_list = Vec::new();
        live_indices.push(*id);
        for (_near_point_distance, near_point_id) in
            kdtree.within(&point, squared_distance, &squared_euclidean)?
        {
            remove_list.push(*near_point_id);
        }
        for key in remove_list {
            points_map.remove(&key);
        }
    }
    Ok(live_indices
        .into_iter()
        .map(|i| points[i].to_vec())
        .collect())
}

#[test]
fn det_4x4_test() {
    let r1 = vec![1., 2., 3., 4.];
    let r2 = vec![5., 6., 7., 8.];
    let r3 = vec![9., 10., 11., 12.];
    let r4 = vec![13., 14., 15., 16.];
    assert_eq!(det(&[r1, r2, r3, r4]), 0.);
    let r1 = vec![1., 2., 3., 1.];
    let r2 = vec![6., 6., 7., 1.];
    let r3 = vec![9., 11., 11., 1.];
    let r4 = vec![15., 14., 15., 1.];
    assert_eq!(det(&[r1, r2, r3, r4]), -4.);
    let r1 = vec![1., 1., 1., 1.];
    let r2 = vec![1., 2., 2., 1.];
    let r3 = vec![1., 2., 3., 1.];
    let r4 = vec![1., 1., 1., 2.];
    assert_eq!(det(&[r1, r2, r3, r4]), 1.);
}
