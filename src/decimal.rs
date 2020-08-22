use super::util::*;
use std::fmt::Display;
use std::collections::{BTreeMap, BTreeSet};
use bigdecimal::{BigDecimal,Zero,One,FromPrimitive};
use std::str::FromStr;
use super::convex::{ConvexHull, Facet, ErrorKind};
use num_bigint::Sign;
use num_traits::Float;

impl Facet<BigDecimal> {
    fn new_bd(points: &[Vec<BigDecimal>], indices: &[usize]) -> Self {
        let points_of_facet: Vec<_> = indices.iter().map(|i| points[*i].to_vec()).collect();
        let normal = facet_normal(&points_of_facet);
        let origin = normal
            .iter()
            .zip(points_of_facet[0].iter())
            .map(|(a, b)| a * b)
            .fold(BigDecimal::zero(), |sum, x| sum + x);
        Self {
            indices: indices.to_vec(),
            outside_points: Vec::new(),
            neighbor_facets: Vec::new(),
            normal,
            origin,
        }
    }
}

impl ConvexHull<BigDecimal> {
    pub fn try_new_bd_from_float<T:Float+Display>(points: &[Vec<T>], max_iter: Option<usize>) -> Result<Self, ErrorKind>{
        let mut bd_points = Vec::new();
        for point in points{
            bd_points.push(float_to_bigdecimal(point));
        }
        Self::try_new_bd(&bd_points, max_iter)
    }
    pub fn try_new_bd(points: &[Vec<BigDecimal>], max_iter: Option<usize>) -> Result<Self, ErrorKind> {
        let num_points = points.len();
        if num_points == 0 {
            return Err(ErrorKind::Empty);
        }
        let dim = points[0].len();
        if dim < 2 {
            return Err(ErrorKind::LessThanTwoDim);
        }
        if !is_same_dimension(&points) {
            return Err(ErrorKind::WrongDimension);
        }
        if num_points <= dim || is_degenerate(&points) {
            return Err(ErrorKind::Degenerated);
        }
        // create simplex of dim+1 points
        let mut c_hull = Self::create_simplex_bd(points)?;
        // main quick hull algorithm
        c_hull.update_bd(max_iter)?;
        // shrink
        c_hull.remove_unused_points();
        if c_hull.points.len() <= dim {
            //degenerate convex hull is generated
            return Err(ErrorKind::Degenerated);
        }
        Ok(c_hull)
    }

    fn create_simplex_bd(points: &[Vec<BigDecimal>]) -> Result<Self, ErrorKind> {
        let indices_set = select_vertices_for_simplex(&points)?;
        let dim = points[0].len();
        let mut facet_add_count = 0;
        let mut facets = BTreeMap::new();
        for i_facet in 0..dim + 1 {
            let mut facet = Vec::new();
            // create facet
            for j in 0..dim + 1 {
                if j == i_facet {
                    continue;
                }
                facet.push(indices_set[j]);
            }
            // check the order of facet's vertices
            let rem_point = indices_set[i_facet];
            let mut mat = Vec::new();
            for point in facet.iter() {
                let mut row = points[*point].to_vec();
                row.push(BigDecimal::one());
                mat.push(row);
            }
            let mut row = points[rem_point].to_vec();
            row.push(BigDecimal::one());
            mat.push(row);
            if det(&mat) < BigDecimal::zero() {
                facet.swap(0, 1);
            }
            debug_assert_eq!(dim, facet.len(), "number of facet's vartices should be dim");
            let facet = Facet::new_bd(points, &facet);
            facets.insert(facet_add_count, facet);
            facet_add_count += 1;
        }
        // link neighbors
        let simplex_facet_key: Vec<_> = facets.keys().map(|k| *k).collect();
        for (key, facet) in &mut facets.iter_mut() {
            facet
                .outside_points
                .sort_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap());
            for neighbors_key in simplex_facet_key
                .iter()
                .filter(|neighbor_key| *neighbor_key != key)
            {
                facet.neighbor_facets.push(*neighbors_key);
            }
        }
        let simplex = Self {
            points: points.to_vec(),
            facets,
        };
        Ok(simplex)
    }

    fn update_bd(&mut self, max_iter: Option<usize>) -> Result<(), ErrorKind> {
        let dim = self.points[0].len();
        let mut facet_add_count = *self.facets.iter().last().map(|(k, _v)| k).unwrap() + 1;
        let mut num_iter = 0;
        let mut assigned_point_indices: BTreeSet<usize> = BTreeSet::new();
        for facet in self.facets.values() {
            for index in &facet.indices {
                assigned_point_indices.insert(*index);
            }
        }
        // initialize outside points
        for (_key, facet) in &mut self.facets.iter_mut() {
            for (i, _point) in self.points.iter().enumerate() {
                if assigned_point_indices.contains(&i) {
                    continue;
                }
                let pos = position_from_facet(&self.points, facet, i);
                //if pos.sign() == Sign::Plus {
                if (pos.clone() - BigDecimal::from_str("0.00001").unwrap()).sign() == Sign::Plus {
                    facet.outside_points.push((i, pos));
                }
            }
        }
        // main algorithm of quick hull
        while let Some((key, facet)) = self
            .facets
            .iter()
            .find(|(_, facet)| !facet.outside_points.is_empty())
            .map(|(a, b)| (*a, b.clone()))
        {
            if let Some(max_iter) = max_iter {
                if num_iter >= max_iter {
                    break;
                }
            }
            num_iter += 1;
            // select the furthest point
            let (furthest_point_index, _) = *facet.outside_points.last().unwrap();
            // initialize visible set
            let visible_set = initialize_visible_set(
                &self.points,
                furthest_point_index,
                &self.facets,
                key,
                &facet,
            );
            // get horizon
            let horizon = get_horizon(&visible_set, &self.facets, dim)?;
            // create new facet
            let mut new_keys = Vec::new();
            for (ridge, unvisible) in horizon {
                let mut new_facet = vec![furthest_point_index];
                assigned_point_indices.insert(furthest_point_index);
                for point in ridge {
                    new_facet.push(point);
                    assigned_point_indices.insert(point);
                }
                if new_facet.len() != dim {
                    return Err(ErrorKind::RoundOffError(
                        "number of new facet's vertices should be dim".to_string(),
                    ));
                }

                let mut new_facet = Facet::new_bd(&self.points, &new_facet);
                new_facet.neighbor_facets.push(unvisible);
                let new_key = facet_add_count;
                facet_add_count += 1;
                self.facets.insert(new_key, new_facet);
                let unvisible_faset = self.facets.get_mut(&unvisible).unwrap();
                unvisible_faset.neighbor_facets.push(new_key);
                new_keys.push(new_key);
            }
            if new_keys.len() < dim {
                return Err(ErrorKind::RoundOffError(
                    "number of new facets should be grater than dim".to_string(),
                ));
            }
            // facet link its neighbor
            for (i, key_a) in new_keys.iter().enumerate() {
                let points_of_new_facet_a: BTreeSet<_> = self
                    .facets
                    .get(key_a)
                    .unwrap()
                    .indices
                    .iter()
                    .map(|k| *k)
                    .collect();
                for key_b in new_keys.iter().skip(i + 1) {
                    let points_of_new_facet_b: BTreeSet<_> = self
                        .facets
                        .get(key_b)
                        .unwrap()
                        .indices
                        .iter()
                        .map(|k| *k)
                        .collect();
                    let num_intersection_points = points_of_new_facet_a
                        .intersection(&points_of_new_facet_b)
                        .collect::<Vec<_>>()
                        .len();
                    if num_intersection_points == dim - 1 {
                        {
                            let facet_a = self.facets.get_mut(key_a).unwrap();
                            facet_a.neighbor_facets.push(*key_b);
                        }
                        let facet_b = self.facets.get_mut(key_b).unwrap();
                        facet_b.neighbor_facets.push(*key_a);
                    }
                }
                let facet_a = self.facets.get(key_a).unwrap();
                if facet_a.neighbor_facets.len() != dim {
                    return Err(ErrorKind::RoundOffError(
                        "number of neighbors should be dim".to_string(),
                    ));
                }
            }
            // check the order of new facet's vertices
            for new_key in &new_keys {
                let new_facet = self.facets.get(new_key).unwrap();
                let mut degenerate = true;
                for assigned_point_index in &assigned_point_indices {
                    let position =
                        position_from_facet(&self.points, &new_facet, *assigned_point_index);
                    //if position.sign() == Sign::NoSign {
                    if (position.clone() - BigDecimal::from_str("0.00001").unwrap()).abs()  <= BigDecimal::from_str("0.00002").unwrap() {
                        continue;
                    } else if position.sign() == Sign::Plus {
                        let new_facet = self.facets.get_mut(new_key).unwrap();
                        new_facet.indices.swap(0, 1);
                        new_facet.normal = new_facet.normal.iter().map(|x| -x).collect();
                        new_facet.origin = -new_facet.origin.clone();
                        degenerate = false;
                        break;
                    } else {
                        degenerate = false;
                        break;
                    }
                }
                if degenerate {
                    //degenerate facet was generated
                    return Err(ErrorKind::Degenerated);
                }
            }
            // set outside points for each new facets
            let mut visible_facets = Vec::new();
            for visible in &visible_set {
                visible_facets.push(self.facets.get(&visible).unwrap().clone());
            }
            for new_key in &new_keys {
                let new_facet = self.facets.get_mut(&new_key).unwrap();
                let mut checked_point_set = BTreeSet::new();
                for visible_facet in &visible_facets {
                    for (outside_point_index, _) in visible_facet.outside_points.iter() {
                        if assigned_point_indices.contains(outside_point_index) {
                            continue;
                        }
                        if checked_point_set.contains(outside_point_index) {
                            continue;
                        } else {
                            checked_point_set.insert(outside_point_index);
                        }
                        let pos =
                            position_from_facet(&self.points, new_facet, *outside_point_index);
                        if (pos.clone() - BigDecimal::from_str("0.00001").unwrap()).sign() == Sign::Plus {
                            new_facet.outside_points.push((*outside_point_index, pos));
                        }
                    }
                }
                new_facet
                    .outside_points
                    .sort_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap());
            }
            // delete the visible facets
            for visible in visible_set {
                let visible_facet = self.facets.get(&visible).unwrap().clone();
                for neighbor_key in visible_facet.neighbor_facets {
                    let neighbor = self.facets.get_mut(&neighbor_key).unwrap();
                    let index = neighbor
                        .neighbor_facets
                        .iter()
                        .enumerate()
                        .find(|(_, k)| **k == visible)
                        .map(|(i, _)| i)
                        .unwrap();
                    neighbor.neighbor_facets.swap_remove(index);
                }
                self.facets.remove(&visible);
            }
        }
        Ok(())
    }

    pub fn add_points_bd(&mut self, points: &[Vec<BigDecimal>]) -> Result<(), ErrorKind> {
        let mut points = points.to_vec();
        self.points.append(&mut points);
        if !is_same_dimension(&self.points) {
            return Err(ErrorKind::WrongDimension);
        }
        self.update_bd(None)?;
        self.remove_unused_points();
        if self.points.len() <= points[0].len() {
            //degenerate convex hull is generated
            return Err(ErrorKind::Degenerated);
        }
        Ok(())
    }
    pub fn volume_bd(&self) -> BigDecimal {
        let dim = self.points[0].len();
        let (c_hull_vertices, c_hull_indices) = self.vertices_indices();
        let mut reference_point = c_hull_vertices[c_hull_indices[0]].to_vec();
        reference_point.push(BigDecimal::one());
        let mut volume = BigDecimal::zero();
        for i in (dim..c_hull_indices.len()).step_by(dim) {
            let mut mat = Vec::new();
            for j in 0..dim {
                let mut row = c_hull_vertices[c_hull_indices[i + j]].to_vec();
                row.push(BigDecimal::one());
                mat.push(row);
            }
            mat.push(reference_point.to_vec());
            volume = volume + det(&mat);
        }
        let factorial = {
            let mut result = 1;
            let mut n = dim;
            while n > 1 {
                result = result * n;
                n = n - 1;
            }
            result
        };
        volume / BigDecimal::from_usize(factorial).unwrap()
    }
    pub fn area_bd(&self) -> BigDecimal {
        let dim = self.points[0].len();
        let (c_hull_vertices, c_hull_indices) = self.vertices_indices();
        let mut total_area = BigDecimal::zero();
        for i in (0..c_hull_indices.len()).step_by(dim) {
            let mut points = Vec::new();
            for j in 0..dim {
                points.push(c_hull_vertices[c_hull_indices[i + j]].to_vec());
            }
            total_area = total_area + facet_area(&points);
        }
        total_area
    }
    pub fn support_point_bd(&self, direction: &[BigDecimal]) -> Result<Vec<BigDecimal>, ErrorKind> {
        if self.points[0].len() != direction.len() {
            return Err(ErrorKind::WrongDimension);
        }
        let (mut facet_key, _) = &self.facets.iter().next().unwrap();
        let facet = &self.facets[facet_key];
        let dot = facet
            .normal
            .iter()
            .zip(direction.iter())
            .map(|(a, b)| a * b)
            .fold(BigDecimal::zero(), |sum, x| sum + x);
        let mut max_dot = dot;
        let mut max_dot_facet_key = facet_key;
        loop {
            let facet = &self.facets[&facet_key];
            for neighbor_key in &facet.neighbor_facets {
                let neighbor = &self.facets[neighbor_key];
                let neighbor_dot = neighbor
                    .normal
                    .iter()
                    .zip(direction.iter())
                    .map(|(a, b)| a * b)
                    .fold(BigDecimal::zero(), |sum, x| sum + x);
                if max_dot < neighbor_dot {
                    max_dot = neighbor_dot;
                    max_dot_facet_key = neighbor_key;
                }
            }
            if facet_key == max_dot_facet_key {
                break;
            } else {
                facet_key = max_dot_facet_key;
            }
        }
        let mut max = BigDecimal::zero();
        let mut max_index = 0;
        for (i,index) in self.facets[&facet_key].indices.iter().enumerate() {
            let dot = self.points[*index]
                .iter()
                .zip(direction.iter())
                .map(|(a, b)| a * b)
                .fold(BigDecimal::zero(), |sum, x| sum + x);
            if i==0 || max < dot {
                max = dot;
                max_index = *index;
            }
        }
        Ok(self.points[max_index].to_vec())
    }
}



fn select_vertices_for_simplex(points: &[Vec<BigDecimal>]) -> Result<Vec<usize>, ErrorKind> {
    // try find the min max point
    let min_max_index_each_axis = min_max_index_each_axis(points);
    let mut vertex_indices_for_simplex = Vec::new();
    for (k, (min_index, max_index)) in min_max_index_each_axis.iter().enumerate() {
        vertex_indices_for_simplex.push(*max_index);
        if k == 0 {
            vertex_indices_for_simplex.push(*min_index);
        }
    }
    // if the point is degenerate, the non-degenerate point is selected again
    if is_degenerate(
        &vertex_indices_for_simplex
            .iter()
            .map(|i| points[*i].to_vec())
            .collect::<Vec<_>>(),
    ) {
        if let Some(indices) = non_degenerate_indices(points) {
            vertex_indices_for_simplex = indices;
        } else {
            return Err(ErrorKind::Degenerated);
        }
    }
    debug_assert_eq!(
        points[0].len() + 1,
        vertex_indices_for_simplex.len(),
        "number of simplex's vertices should be dim+1"
    );
    Ok(vertex_indices_for_simplex)
}

// get visible facet viewed from furthest point
fn initialize_visible_set(
    points: &[Vec<BigDecimal>],
    furthest_point_index: usize,
    facets: &BTreeMap<usize, Facet<BigDecimal>>,
    faset_key: usize,
    facet: &Facet<BigDecimal>,
) -> BTreeSet<usize> {
    let mut visible_set = BTreeSet::new();
    visible_set.insert(faset_key);
    let mut neighbor_stack: Vec<_> = facet.neighbor_facets.iter().map(|k| *k).collect();
    let mut visited_neighbor = BTreeSet::new();
    while let Some(neighbor_key) = neighbor_stack.pop() {
        if visited_neighbor.contains(&neighbor_key) {
            continue;
        } else {
            visited_neighbor.insert(neighbor_key);
        }
        let neighbor = facets.get(&neighbor_key).unwrap();
        let pos = position_from_facet(points, neighbor, furthest_point_index);
        if (pos.clone() - BigDecimal::from_str("0.00001").unwrap()).sign() == Sign::Plus {
            visible_set.insert(neighbor_key);
            neighbor_stack.append(&mut neighbor.neighbor_facets.iter().map(|k| *k).collect());
        }
    }
    visible_set
}

fn get_horizon(
    visible_set: &BTreeSet<usize>,
    facets: &BTreeMap<usize, Facet<BigDecimal>>,
    dim: usize, // assertion use only
) -> Result<Vec<(Vec<usize>, usize)>, ErrorKind> {
    let mut horizon = Vec::new();
    for visible_key in visible_set {
        let visible_facet = facets.get(visible_key).unwrap();
        let points_of_visible_facet: BTreeSet<_> =
            visible_facet.indices.iter().map(|i| *i).collect();
        if dim != points_of_visible_facet.len() {
            return Err(ErrorKind::RoundOffError(
                "number of visible facet's vartices should be dim".to_string(),
            ));
        }

        for neighbor_key in &visible_facet.neighbor_facets {
            // if neighbor is unvisible
            if !visible_set.contains(neighbor_key) {
                let unvisible_neighbor = facets.get(neighbor_key).unwrap();
                let points_of_unvisible_neighbor: BTreeSet<_> =
                    unvisible_neighbor.indices.iter().map(|i| *i).collect();
                if dim != points_of_unvisible_neighbor.len() {
                    return Err(ErrorKind::RoundOffError(
                        "number of unvisible facet's vartices should be dim".to_string(),
                    ));
                }

                let horizon_ridge: Vec<_> = points_of_unvisible_neighbor
                    .intersection(&points_of_visible_facet)
                    .map(|key| *key)
                    .collect();
                if dim - 1 != horizon_ridge.len() {
                    return Err(ErrorKind::RoundOffError(
                        "number of ridge's vartices should be dim-1".to_string(),
                    ));
                }
                horizon.push((horizon_ridge, *neighbor_key));
            }
        }
    }
    if horizon.len() < dim{
        return Err(ErrorKind::RoundOffError("horizon len < dim".to_string()));
    }
    Ok(horizon)
}

fn position_from_facet(points: &[Vec<BigDecimal>], facet: &Facet<BigDecimal>, point_index: usize) -> BigDecimal {
    let point = points[point_index].to_vec();
    let origin = facet.origin.clone();
    let pos = facet
        .normal
        .iter()
        .zip(point.iter())
        .map(|(a, b)| a * b)
        .fold(BigDecimal::zero(), |sum, x| sum + x);
    pos - origin
}

fn is_degenerate(points: &[Vec<BigDecimal>]) -> bool {
    let dim = points[0].len();
    let ex_vec: Vec<Vec<_>> = points
        .iter()
        .map(|v| {
            let mut v = v.to_vec();
            v.push(BigDecimal::one());
            v
        })
        .collect();
    let num = ex_vec.len();
    if dim >= num {
        return true;
    }
    let mut mat = Vec::new();
    for i in 0..dim + 1 {
        let mut row = Vec::new();
        for j in 0..dim + 1 {
            let mut c = BigDecimal::zero();
            for k in 0..num {
                c = c + ex_vec[k][i].clone() * ex_vec[k][j].clone();
            }
            row.push(c);
        }
        mat.push(row);
    }
    //if det(&mat).sign() == Sign::NoSign {
    if (det(&mat)-BigDecimal::from_str("0.00001").unwrap()).abs() <= BigDecimal::from_str("0.00002").unwrap() {
        true
    } else {
        false
    }
}

fn non_degenerate_indices(vertices: &[Vec<BigDecimal>]) -> Option<Vec<usize>> {
    // look for points that are not degenerate for simplex using the Gram-Schmidt method
    let dim = vertices[0].len();
    if dim >= vertices.len() {
        return None;
    }
    let mut indices = vec![0];
    let mut first_axis = vertices[1]
        .iter()
        .zip(vertices[0].iter())
        .map(|(a, b)| a - b)
        .collect::<Vec<_>>();
    let mut sq_norm = first_axis.iter().fold(BigDecimal::zero(), |sum, x| sum + x.square());
    let mut j = 1;
    //while (sq_norm).sign() == Sign::NoSign {
    while (sq_norm.clone() -BigDecimal::from_str("0.00001").unwrap()).abs() <= BigDecimal::from_str("0.00002").unwrap() {
        j += 1;
        first_axis = vertices[j]
            .iter()
            .zip(vertices[0].iter())
            .map(|(a, b)| a - b)
            .collect();
        sq_norm = first_axis.iter().fold(BigDecimal::zero(), |sum, x| sum + x.square());
    }
    let norm = sq_norm.sqrt().unwrap();
    let first_axis: Vec<_> = first_axis.into_iter().map(|x| x / norm.clone()).collect();
    let mut axes = vec![first_axis];
    indices.push(j);
    for i in j + 1..vertices.len() {
        let vector = vertices[i]
            .iter()
            .zip(vertices[0].iter())
            .map(|(a, b)| a - b);
        let norm = vector.clone().fold(BigDecimal::zero(), |sum, x| sum + x.square()).sqrt().unwrap();
        if (norm.clone() -BigDecimal::from_str("0.00001").unwrap()).abs() <= BigDecimal::from_str("0.00002").unwrap() {
            continue;
        }
        let unit_vector: Vec<_> = vector.map(|x| x / norm.clone()).collect();
        let mut rem = unit_vector.to_vec();
        for axis in &axes {
            let coef = axis
                .iter()
                .zip(unit_vector.iter())
                .map(|(a, b)| a * b)
                .fold(BigDecimal::zero(), |sum, x| sum + x);
            rem = rem
                .iter()
                .zip(axis.iter())
                .map(|(rem, axis)| rem - coef.clone() * axis)
                .collect();
        }

        let rem_sq_norm = rem.iter().fold(BigDecimal::zero(), |sum, x| sum + x.square());
        if (rem_sq_norm.clone()-BigDecimal::from_str("0.00001").unwrap()).sign() == Sign::Minus {
            continue;
        }
        let rem_norm = rem_sq_norm.sqrt().unwrap();
        let new_axis: Vec<_> = rem.into_iter().map(|x| x / rem_norm.clone()).collect();
        axes.push(new_axis);
        indices.push(i);
        if axes.len() == dim {
            return Some(indices);
        }
    }
    None
}

fn facet_area(points_of_facet: &[Vec<BigDecimal>]) -> BigDecimal {
    let num_points = points_of_facet.len();
    let dim = points_of_facet[0].len();
    debug_assert_eq!(num_points, dim);
    let mut mat = Vec::new();
    for i in 1..num_points {
        let row: Vec<_> = points_of_facet[i]
            .iter()
            .zip(points_of_facet[i - 1].iter())
            .map(|(a, b)| a - b)
            .collect();
        mat.push(row);
    }
    mat.push(facet_normal(points_of_facet));
    let factorial = {
        let mut result = 1;
        let mut n = dim - 1;
        while n > 1 {
            result = result * n;
            n = n - 1;
        }
        result
    };
    det(&mat).abs() / BigDecimal::from_usize(factorial).unwrap()
}

fn facet_normal(points_of_facet: &[Vec<BigDecimal>]) -> Vec<BigDecimal> {
    let num_points = points_of_facet.len();
    let dim = points_of_facet[0].len();
    debug_assert_eq!(num_points, dim);
    let mut vectors = Vec::new();
    for i in 1..num_points {
        let vector: Vec<_> = points_of_facet[i]
            .iter()
            .zip(points_of_facet[i - 1].iter())
            .map(|(a, b)| a - b)
            .collect();
        vectors.push(vector);
    }
    let mut sign = BigDecimal::one();
    let mut normal = Vec::new();
    for i in 0..dim {
        let mut mat = Vec::new();
        for vector in vectors.iter() {
            let mut column = Vec::new();
            for (j, element) in vector.iter().enumerate() {
                if j == i {
                    continue;
                }
                column.push(element.clone());
            }
            mat.push(column);
        }
        let cofactor = det(&mat);
        normal.push(sign.clone() * cofactor);
        sign = -sign;
    }
    let norm = normal.iter().fold(BigDecimal::zero(), |a, x| a + x.square()).sqrt().unwrap();
    let normalized_normal = normal.into_iter().map(|x| x / norm.clone()).collect();
    normalized_normal
}

fn float_to_bigdecimal<T: Float + Display>(vec: &[T]) -> Vec<BigDecimal>{
    let mut decimal_vec = Vec::new();
    for v in vec{
        decimal_vec.push(BigDecimal::from_str(&format!("{:.10}", v)).unwrap());
    }
    decimal_vec
}