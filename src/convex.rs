use kdtree::distance::squared_euclidean;
use kdtree::KdTree;
use num_traits::Float;
use std::collections::{BTreeMap, BTreeSet};
use std::error::Error;
use std::fmt;

#[derive(Debug, Clone)]
struct Facet<T: Float> {
    indices: Vec<usize>,
    outside_points: Vec<(usize, T)>,
    neighbor_facets: Vec<usize>,
    normal: Vec<T>,
    origin: T,
}

impl<T: Float> Facet<T> {
    fn new(points: &[Vec<T>], indices: &[usize]) -> Self {
        let points_of_facet: Vec<_> = indices.iter().map(|i| points[*i].to_vec()).collect();
        let normal = facet_normal(&points_of_facet);
        let origin = normal
            .iter()
            .zip(points_of_facet[0].iter())
            .map(|(a, b)| *a * *b)
            .fold(T::zero(), |sum, x| sum + x);
        Self {
            indices: indices.to_vec(),
            outside_points: Vec::new(),
            neighbor_facets: Vec::new(),
            normal,
            origin,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum ErrorKind {
    Empty,
    LessThanTwoDim,
    Degenerated,
    WrongDimension,
    RoundOffError(String),
}

impl fmt::Display for ErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self {
            ErrorKind::Empty => write!(f, "empty"),
            ErrorKind::LessThanTwoDim => write!(f, "less than two dimention"),
            ErrorKind::Degenerated => write!(f, "degenerated"),
            ErrorKind::WrongDimension => write!(f, "wrong dimension"),
            ErrorKind::RoundOffError(msg) => {
                write!(f, "erroneous results by roundoff error: {}", msg)
            }
        }
    }
}
impl Error for ErrorKind {}

#[derive(Clone, Debug)]
pub struct ConvexHull<T: Float> {
    points: Vec<Vec<T>>,
    facets: BTreeMap<usize, Facet<T>>,
}

impl<T: Float> ConvexHull<T> {
    pub fn try_new(points: &[Vec<T>], max_iter: Option<usize>) -> Result<Self, ErrorKind> {
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
        // remove nearby points
        let points = &remove_nearby_points(&points, T::epsilon() * T::from(15).unwrap());
        if num_points <= dim || is_degenerate(&points) {
            return Err(ErrorKind::Degenerated);
        }
        // create simplex of dim+1 points
        let mut c_hull = Self::create_simplex(points)?;
        // main quick hull algorithm
        c_hull.update(max_iter)?;
        // shrink
        c_hull.remove_unused_points();
        if c_hull.points.len() <= dim {
            //degenerate convex hull is generated
            return Err(ErrorKind::Degenerated);
        }
        Ok(c_hull)
    }

    fn create_simplex(points: &[Vec<T>]) -> Result<Self, ErrorKind> {
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
                row.push(T::one());
                mat.push(row);
            }
            let mut row = points[rem_point].to_vec();
            row.push(T::one());
            mat.push(row);
            if det(&mat) < T::zero() {
                facet.swap(0, 1);
            }
            debug_assert_eq!(dim, facet.len(), "number of facet's vartices should be dim");
            let facet = Facet::new(points, &facet);
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

    fn update(&mut self, max_iter: Option<usize>) -> Result<(), ErrorKind> {
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
                if pos > T::epsilon() * T::from(200).unwrap() {
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

                let mut new_facet = Facet::new(&self.points, &new_facet);
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
                let mut mat = Vec::new();
                for index in &new_facet.indices{
                    let mut row = self.points[*index].to_vec();
                    row.push(T::one());
                    mat.push(row);
                }
                for assigned_point_index in &assigned_point_indices {
                    let mut row = self.points[*assigned_point_index].to_vec();
                    row.push(T::one());
                    mat.push(row);
                    let det = det(&mat);
                    if  det.abs() == T::zero() {
                        mat.pop();
                        continue;
                    } else if det < T::zero() {
                        let new_facet = self.facets.get_mut(new_key).unwrap();
                        new_facet.indices.swap(0, 1);
                        new_facet.normal = new_facet.normal.iter().map(|x| -(*x)).collect();
                        new_facet.origin = -new_facet.origin;
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
                        if pos > T::epsilon() * T::from(200).unwrap() {
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
        if self.volume() < T::epsilon(){
            Err(ErrorKind::RoundOffError("volume zero".to_string()))
        }else{
            Ok(())
        }
    }

    pub fn add_points(&mut self, points: &[Vec<T>]) -> Result<(), ErrorKind> {
        let mut points = points.to_vec();
        self.points.append(&mut points);
        if !is_same_dimension(&self.points) {
            return Err(ErrorKind::WrongDimension);
        }
        self.update(None)?;
        self.remove_unused_points();
        if self.points.len() <= points[0].len() {
            //degenerate convex hull is generated
            return Err(ErrorKind::Degenerated);
        }
        Ok(())
    }
    pub fn vertices_indices(&self) -> (Vec<Vec<T>>, Vec<usize>) {
        let mut indices = Vec::new();
        for facet in self.facets.values() {
            for i in &facet.indices {
                indices.push(*i);
            }
        }
        (self.points.to_vec(), indices)
    }
    fn remove_unused_points(&mut self) {
        let mut indices_list = BTreeSet::new();
        for facet in self.facets.values() {
            for i in &facet.indices {
                indices_list.insert(*i);
            }
        }
        let indices_list: BTreeMap<usize, usize> = indices_list
            .into_iter()
            .enumerate()
            .map(|(i, index)| (index, i))
            .collect();
        for facet in self.facets.values_mut() {
            let mut new_facet_indices = Vec::new();
            for i in &facet.indices {
                new_facet_indices.push(*indices_list.get(&i).unwrap());
            }
            std::mem::swap(&mut facet.indices, &mut new_facet_indices);
        }
        let mut vertices = Vec::new();
        for (index, _i) in indices_list.iter() {
            let point = self.points[*index].to_vec();
            vertices.push(point);
        }
        self.points = vertices;
    }
    pub fn volume(&self) -> T {
        let dim = self.points[0].len();
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
            let mut result = 1;
            let mut n = dim;
            while n > 1 {
                result = result * n;
                n = n - 1;
            }
            result
        };
        volume / T::from(factorial).unwrap()
    }
    pub fn area(&self) -> T {
        let dim = self.points[0].len();
        let (c_hull_vertices, c_hull_indices) = self.vertices_indices();
        let mut total_area = T::zero();
        for i in (0..c_hull_indices.len()).step_by(dim) {
            let mut points = Vec::new();
            for j in 0..dim {
                points.push(c_hull_vertices[c_hull_indices[i + j]].to_vec());
            }
            total_area = total_area + facet_area(&points);
        }
        total_area
    }
    pub fn support_point(&self, direction: &[T]) -> Result<Vec<T>, ErrorKind> {
        if self.points[0].len() != direction.len() {
            return Err(ErrorKind::WrongDimension);
        }
        let (mut facet_key, _) = &self.facets.iter().next().unwrap();
        let facet = &self.facets[facet_key];
        let dot = facet
            .normal
            .iter()
            .zip(direction.iter())
            .map(|(a, b)| *a * *b)
            .fold(T::zero(), |sum, x| sum + x);
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
                    .map(|(a, b)| *a * *b)
                    .fold(T::zero(), |sum, x| sum + x);
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
        let mut max = T::neg_infinity();
        let mut max_index = 0;
        for index in &self.facets[&facet_key].indices {
            let dot = self.points[*index]
                .iter()
                .zip(direction.iter())
                .map(|(a, b)| *a * *b)
                .fold(T::zero(), |sum, x| sum + x);
            if max < dot {
                max = dot;
                max_index = *index;
            }
        }
        Ok(self.points[max_index].to_vec())
    }
}

fn min_max_index_each_axis<T: Float>(points: &[Vec<T>]) -> Vec<(usize, usize)> {
    let dim = points[0].len();
    let mut min_index = vec![0; dim];
    let mut max_index = vec![0; dim];
    let mut min = vec![Float::infinity(); dim];
    let mut max = vec![Float::neg_infinity(); dim];
    for (index, point) in points.iter().enumerate() {
        for (j, element) in point.iter().enumerate() {
            if *element < min[j] {
                min[j] = *element;
                min_index[j] = index;
            }
            if *element > max[j] {
                max[j] = *element;
                max_index[j] = index;
            }
        }
    }
    min_index.into_iter().zip(max_index.into_iter()).collect()
}

fn select_vertices_for_simplex<T: Float>(points: &[Vec<T>]) -> Result<Vec<usize>, ErrorKind> {
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
fn initialize_visible_set<T: Float>(
    points: &[Vec<T>],
    furthest_point_index: usize,
    facets: &BTreeMap<usize, Facet<T>>,
    faset_key: usize,
    facet: &Facet<T>,
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
        if pos > T::epsilon() * T::from(200).unwrap() {
            visible_set.insert(neighbor_key);
            neighbor_stack.append(&mut neighbor.neighbor_facets.iter().map(|k| *k).collect());
        }
    }
    visible_set
}

fn get_horizon<T: Float>(
    visible_set: &BTreeSet<usize>,
    facets: &BTreeMap<usize, Facet<T>>,
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

fn position_from_facet<T: Float>(points: &[Vec<T>], facet: &Facet<T>, point_index: usize) -> T {
    let point = points[point_index].to_vec();
    let origin = facet.origin;
    let pos = facet
        .normal
        .iter()
        .zip(point.iter())
        .map(|(a, b)| *a * *b)
        .fold(T::zero(), |sum, x| sum + x);
    pos - origin
}

fn is_same_dimension<T>(points: &[Vec<T>]) -> bool {
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

fn is_degenerate<T: Float>(points: &[Vec<T>]) -> bool {
    let dim = points[0].len();
    let ex_vec: Vec<Vec<_>> = points
        .iter()
        .map(|v| {
            let mut v = v.to_vec();
            v.push(T::one());
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
            let mut c = T::zero();
            for k in 0..num {
                c = c + ex_vec[k][i] * ex_vec[k][j];
            }
            row.push(c);
        }
        mat.push(row);
    }
    if det(&mat).abs() <= T::epsilon() * T::from(200).unwrap() {
        true
    } else {
        false
    }
}

fn non_degenerate_indices<T: Float>(vertices: &[Vec<T>]) -> Option<Vec<usize>> {
    // look for points that are not degenerate for simplex using the Gram-Schmidt method
    let dim = vertices[0].len();
    if dim >= vertices.len() {
        return None;
    }
    let mut indices = vec![0];
    let mut first_axis = vertices[1]
        .iter()
        .zip(vertices[0].iter())
        .map(|(a, b)| *a - *b)
        .collect::<Vec<_>>();
    let mut sq_norm = first_axis.iter().fold(T::zero(), |sum, x| sum + *x * *x);
    let mut j = 1;
    while sq_norm <= T::epsilon() * T::from(200).unwrap() {
        j += 1;
        first_axis = vertices[j]
            .iter()
            .zip(vertices[0].iter())
            .map(|(a, b)| *a - *b)
            .collect();
        sq_norm = first_axis.iter().fold(T::zero(), |sum, x| sum + *x * *x);
    }
    let norm = sq_norm.sqrt();
    let first_axis: Vec<_> = first_axis.into_iter().map(|x| x / norm).collect();
    let mut axes = vec![first_axis];
    indices.push(j);
    for i in j + 1..vertices.len() {
        let vector = vertices[i]
            .iter()
            .zip(vertices[0].iter())
            .map(|(a, b)| *a - *b);
        let norm = vector.clone().fold(T::zero(), |sum, x| sum + x * x).sqrt();
        let unit_vector: Vec<_> = vector.map(|x| x / norm).collect();
        let mut rem = unit_vector.to_vec();
        for axis in &axes {
            let coef = axis
                .iter()
                .zip(unit_vector.iter())
                .map(|(a, b)| *a * *b)
                .fold(T::zero(), |sum, x| sum + x);
            rem = rem
                .iter()
                .zip(axis.iter())
                .map(|(rem, axis)| *rem - coef * *axis)
                .collect();
        }

        let rem_sq_norm = rem.iter().fold(T::zero(), |sum, x| sum + *x * *x);
        if rem_sq_norm <= T::epsilon() * T::from(200).unwrap() {
            continue;
        }
        let rem_norm = rem_sq_norm.sqrt();
        let new_axis: Vec<_> = rem.into_iter().map(|x| x / rem_norm).collect();
        axes.push(new_axis);
        indices.push(i);
        if axes.len() == dim {
            return Some(indices);
        }
    }
    None
}

fn remove_nearby_points<T: Float>(points: &[Vec<T>], squared_distance: T) -> Vec<Vec<T>> {
    let mut points_map = BTreeMap::new();
    let dim = points[0].len();
    let mut kdtree = KdTree::new(dim);
    for (i, point) in points.iter().enumerate() {
        points_map.insert(i, point);
        kdtree.add(point, i).unwrap();
    }
    let mut live_indices = Vec::new();
    while let Some((id, point)) = points_map.iter().next() {
        let mut remove_list = Vec::new();
        live_indices.push(*id);
        for (_near_point_distance, near_point_id) in kdtree
            .within(&point, squared_distance, &squared_euclidean)
            .unwrap()
        {
            remove_list.push(*near_point_id);
        }
        for key in remove_list {
            points_map.remove(&key);
        }
    }
    live_indices
        .into_iter()
        .map(|i| points[i].to_vec())
        .collect()
}

fn facet_area<T: Float>(points_of_facet: &[Vec<T>]) -> T {
    let num_points = points_of_facet.len();
    let dim = points_of_facet[0].len();
    debug_assert_eq!(num_points, dim);
    let mut mat = Vec::new();
    for i in 1..num_points {
        let row: Vec<_> = points_of_facet[i]
            .iter()
            .zip(points_of_facet[i - 1].iter())
            .map(|(a, b)| *a - *b)
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
    det(&mat).abs() / T::from(factorial).unwrap()
}

fn facet_normal<T: Float>(points_of_facet: &[Vec<T>]) -> Vec<T> {
    let num_points = points_of_facet.len();
    let dim = points_of_facet[0].len();
    debug_assert_eq!(num_points, dim);
    let mut vectors = Vec::new();
    for i in 1..num_points {
        let vector: Vec<_> = points_of_facet[i]
            .iter()
            .zip(points_of_facet[i - 1].iter())
            .map(|(a, b)| *a - *b)
            .collect();
        vectors.push(vector);
    }
    let mut sign = T::one();
    let mut normal = Vec::new();
    for i in 0..dim {
        let mut mat = Vec::new();
        for vector in vectors.iter() {
            let mut column = Vec::new();
            for (j, element) in vector.iter().enumerate() {
                if j == i {
                    continue;
                }
                column.push(*element);
            }
            mat.push(column);
        }
        let cofactor = det(&mat);
        normal.push(sign * cofactor);
        sign = -sign;
    }
    let norm = normal.iter().fold(T::zero(), |a, &x| a + x * x).sqrt();
    let normalized_normal = normal.into_iter().map(|x| x / norm).collect();
    normalized_normal
}

fn det<T: Float>(m: &[Vec<T>]) -> T {
    let row_dim = m.len();
    let column_dim = m[0].len();
    assert_eq!(row_dim, column_dim);
    match column_dim {
        1 => m[0][0],
        2 => det_2x2(m),
        3 => det_3x3(m),
        4 => det_4x4(m),
        _ => {
            unimplemented!("matrix size should be less than 4 dim to calcurate determinant for now")
        }
    }
}

fn det_2x2<T: Float>(m: &[Vec<T>]) -> T {
    m[0][0] * m[1][1] - m[0][1] * m[1][0]
}

#[rustfmt::skip]
fn det_3x3<T: Float>(m: &[Vec<T>]) -> T {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
  - m[1][0] * (m[0][1] * m[2][2] - m[0][2] * m[2][1])
  + m[2][0] * (m[0][1] * m[1][2] - m[0][2] * m[1][1])
}

#[rustfmt::skip]
fn det_4x4<T: Float>(m: &[Vec<T>]) -> T {
    m[0][0]
        * (   m[1][1] * m[2][2] * m[3][3]
            + m[1][2] * m[2][3] * m[3][1]
            + m[1][3] * m[2][1] * m[3][2]
            - m[1][3] * m[2][2] * m[3][1]
            - m[1][2] * m[2][1] * m[3][3]
            - m[1][1] * m[2][3] * m[3][2])
  - m[1][0]
        * (   m[0][1] * m[2][2] * m[3][3]
            + m[0][2] * m[2][3] * m[3][1]
            + m[0][3] * m[2][1] * m[3][2]
            - m[0][3] * m[2][2] * m[3][1]
            - m[0][2] * m[2][1] * m[3][3]
            - m[0][1] * m[2][3] * m[3][2])
  + m[2][0]
        * (   m[0][1] * m[1][2] * m[3][3]
            + m[0][2] * m[1][3] * m[3][1]
            + m[0][3] * m[1][1] * m[3][2]
            - m[0][3] * m[1][2] * m[3][1]
            - m[0][2] * m[1][1] * m[3][3]
            - m[0][1] * m[1][3] * m[3][2])
  - m[3][0]
        * (   m[0][1] * m[1][2] * m[2][3]
            + m[0][2] * m[1][3] * m[2][1]
            + m[0][3] * m[1][1] * m[2][2]
            - m[0][3] * m[1][2] * m[2][1]
            - m[0][2] * m[1][1] * m[2][3]
            - m[0][1] * m[1][3] * m[2][2])
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

#[test]
fn facet_normal_test() {
    let p1 = vec![-1.0, 0.0, 0.0];
    let p2 = vec![1.0, 0.0, 0.0];
    let p3 = vec![0.0, 1.0, 0.0];
    let normal_z = facet_normal(&[p1, p2, p3]);
    assert_eq!(normal_z, vec!(0.0, 0.0, 1.0));

    let p1 = vec![0.0, -1.0, 0.0];
    let p2 = vec![0.0, 1.0, 0.0];
    let p3 = vec![0.0, 0.0, 1.0];
    let normal_x = facet_normal(&[p1, p2, p3]);
    assert_eq!(normal_x, vec!(1.0, 0.0, 0.0));

    let p1 = vec![0.0, 0.0, -1.0];
    let p2 = vec![0.0, 0.0, 1.0];
    let p3 = vec![1.0, 0.0, 0.0];
    let normal_y = facet_normal(&[p1, p2, p3]);
    assert_eq!(normal_y, vec!(0.0, 1.0, 0.0));
}

#[test]
fn inner_outer_test() {
    let p1 = vec![1.0, 0.0, 0.0];
    let p2 = vec![0.0, 1.0, 0.0];
    let p3 = vec![0.0, 0.0, 1.0];
    let outer_point = vec![0.0, 0.0, 10.0];
    let inner_point = vec![0.0, 0.0, 0.0];
    let whithin_point = vec![1.0, 0.0, 0.0];
    let points = vec![p1, p2, p3, outer_point, inner_point, whithin_point];
    let facet: Facet<f64> = Facet::new(&points, &[0, 1, 2]);
    let outer = position_from_facet(&points, &facet, 3);
    assert!(outer > 0.0);
    let inner = position_from_facet(&points, &facet, 4);
    assert!(inner < 0.0);
    let within = position_from_facet(&points, &facet, 5);
    assert!(within == 0.0);
}

#[test]
fn rectangle_test() {
    let p1 = vec![2.0, 2.0];
    let p2 = vec![4.0, 2.0];
    let p3 = vec![4.0, 4.0];
    let p4 = vec![2.0, 4.0];
    let p5 = vec![3.0, 3.0];

    let rect = ConvexHull::try_new(&[p1, p2, p3, p4, p5], None).unwrap();
    assert_eq!(rect.area(), 8.0);
    assert_eq!(rect.volume(), 4.0);

    let (_v, i) = rect.vertices_indices();
    assert_eq!(i.len(), 4 * 2);
}

#[test]
fn octahedron_test() {
    let p1 = vec![1.0, 0.0, 0.0];
    let p2 = vec![0.0, 1.0, 0.0];
    let p3 = vec![0.0, 0.0, 1.0];
    let p4 = vec![-1.0, 0.0, 0.0];
    let p5 = vec![0.0, -1.0, 0.0];
    let p6 = vec![0.0, 0.0, -1.0];
    let (_v, i) = ConvexHull::try_new(&[p1, p2, p3, p4, p5, p6], None)
        .unwrap()
        .vertices_indices();
    assert_eq!(i.len(), 8 * 3);
}

#[test]
fn octahedron_translation_test() {
    fn translate(points: &Vec<f64>) -> Vec<f64> {
        let d = vec![10.0, 10.0, 10.0];
        points.iter().zip(d.iter()).map(|(a, b)| *a + *b).collect()
    }
    let p1 = vec![1.0, 0.0, 0.0];
    let p2 = vec![0.0, 1.0, 0.0];
    let p3 = vec![0.0, 0.0, 1.0];
    let p4 = vec![-1.0, 0.0, 0.0];
    let p5 = vec![0.0, -1.0, 0.0];
    let p6 = vec![0.0, 0.0, -1.0];
    let points: Vec<_> = [p1, p2, p3, p4, p5, p6]
        .iter()
        .map(|p| translate(p))
        .collect();
    let (_v, i) = ConvexHull::try_new(&points, None)
        .unwrap()
        .vertices_indices();
    assert_eq!(i.len(), 8 * 3);
}

#[test]
fn cube_test() {
    let p1 = vec![1.0, 1.0, 1.0];
    let p2 = vec![1.0, 1.0, -1.0];
    let p3 = vec![1.0, -1.0, 1.0];
    let p4 = vec![1.0, -1.0, -1.0];
    let p5 = vec![-1.0, 1.0, 1.0];
    let p6 = vec![-1.0, 1.0, -1.0];
    let p7 = vec![-1.0, -1.0, 1.0];
    let p8 = vec![-1.0, -1.0, -1.0];
    let (_v, i) = ConvexHull::try_new(&[p1, p2, p3, p4, p5, p6, p7, p8], None)
        .unwrap()
        .vertices_indices();
    assert_eq!(i.len(), 6 * 2 * 3);
}

#[test]
fn cube_area_test() {
    let p1 = vec![2.0, 2.0, 2.0];
    let p2 = vec![2.0, 2.0, 0.0];
    let p3 = vec![2.0, 0.0, 2.0];
    let p4 = vec![2.0, 0.0, 0.0];
    let p5 = vec![0.0, 2.0, 2.0];
    let p6 = vec![0.0, 2.0, 0.0];
    let p7 = vec![0.0, 0.0, 2.0];
    let p8 = vec![0.0, 0.0, 0.0];
    let cube = ConvexHull::try_new(&[p1, p2, p3, p4, p5, p6, p7, p8], None).unwrap();
    assert_eq!(cube.area(), 24.0);
}

#[test]
fn cube_volume_test() {
    let p1 = vec![2.0, 2.0, 2.0];
    let p2 = vec![2.0, 2.0, 0.0];
    let p3 = vec![2.0, 0.0, 2.0];
    let p4 = vec![2.0, 0.0, 0.0];
    let p5 = vec![0.0, 2.0, 2.0];
    let p6 = vec![0.0, 2.0, 0.0];
    let p7 = vec![0.0, 0.0, 2.0];
    let p8 = vec![0.0, 0.0, 0.0];
    let cube = ConvexHull::try_new(&[p1, p2, p3, p4, p5, p6, p7, p8], None).unwrap();
    assert_eq!(cube.volume(), 8.0);
}

#[test]
fn cube_support_point_test() {
    let p1 = vec![1.0, 1.0, 1.0];
    let p2 = vec![1.0, 1.0, 0.0];
    let p3 = vec![1.0, 0.0, 1.0];
    let p4 = vec![1.0, 0.0, 0.0];
    let p5 = vec![0.0, 1.0, 1.0];
    let p6 = vec![0.0, 1.0, 0.0];
    let p7 = vec![0.0, 0.0, 1.0];
    let p8 = vec![0.0, 0.0, 0.0];
    let cube = ConvexHull::try_new(&[p1.to_vec(), p2, p3, p4, p5, p6, p7, p8], None).unwrap();
    assert_eq!(cube.support_point(&vec![0.5, 0.5, 0.5]).unwrap(), p1);
}

#[test]
#[should_panic(expected = "Degenerated")]
fn flat_test() {
    let p1 = vec![1.0, 1.0, 10.0];
    let p2 = vec![1.0, 1.0, 10.0];
    let p3 = vec![1.0, -1.0, 10.0];
    let p4 = vec![1.0, -1.0, 10.0];
    let p5 = vec![-1.0, 1.0, 10.0];
    let p6 = vec![-1.0, 1.0, 10.0];
    let p7 = vec![-1.0, -1.0, 10.0];
    let p8 = vec![-1.0, -1.0, 10.0];
    let (_v, _i) = ConvexHull::try_new(&[p1, p2, p3, p4, p5, p6, p7, p8], None)
        .unwrap()
        .vertices_indices();
}

#[test]
fn simplex_may_degenerate_test() {
    let points = vec![
        vec![1.0, 0.0, 1.0],
        vec![1.0, 1.0, 1.0],
        vec![2.0, 1.0, 0.0],
        vec![2.0, 1.0, 1.0],
        vec![2.0, 0.0, 1.0],
        vec![2.0, 0.0, 0.0],
        vec![1.0, 1.0, 2.0],
        vec![0.0, 1.0, 2.0],
        vec![0.0, 0.0, 2.0],
        vec![1.0, 0.0, 2.0],
    ];
    let (_v, _i) = ConvexHull::try_new(&points, None)
        .unwrap()
        .vertices_indices();
}

#[test]
fn simplex_may_degenerate_test_2() {
    let vertices = vec![
        vec![0., 0., 0.],
        vec![1., 0., 0.],
        vec![1., 0., 1.],
        vec![0., 0., 1.],
        vec![0., 1., 0.],
        vec![1., 1., 0.],
        vec![1., 1., 1.],
        vec![0., 1., 1.],
        vec![2., 1., 0.],
        vec![2., 1., 1.],
        vec![2., 0., 1.],
        vec![2., 0., 0.],
        vec![1., 1., 2.],
        vec![0., 1., 2.],
        vec![0., 0., 2.],
        vec![1., 0., 2.],
    ];
    let indices = vec![4, 5, 1, 11, 1, 5, 1, 11, 10, 10, 2, 1, 5, 8, 11];
    let points = indices
        .iter()
        .map(|i| vertices[*i].to_vec())
        .collect::<Vec<_>>();
    let (_v, _i) = ConvexHull::try_new(&points, None)
        .unwrap()
        .vertices_indices();
}

#[test]
fn sphere_test() {
    fn rot_z(point: &Vec<f64>, angle: f64) -> Vec<f64> {
        let e1 = angle.cos() * point[0] - angle.sin() * point[1];
        let e2 = angle.sin() * point[0] + angle.cos() * point[1];
        let e3 = point[2];
        vec![e1, e2, e3]
    }
    fn rot_x(point: &Vec<f64>, angle: f64) -> Vec<f64> {
        let e1 = point[0];
        let e2 = angle.cos() * point[1] - angle.sin() * point[2];
        let e3 = angle.sin() * point[1] + angle.cos() * point[2];
        vec![e1, e2, e3]
    }
    let mut points = Vec::new();
    let dev = 10;
    let unit_y = vec![0.0, 1.0, 0.0];
    for step_x in 0..dev {
        let angle_x = 2.0 * std::f64::consts::PI * (step_x as f64 / dev as f64);
        let p = rot_x(&unit_y, angle_x);
        for step_z in 0..dev {
            let angle_z = 2.0 * std::f64::consts::PI * (step_z as f64 / dev as f64);
            let p = rot_z(&p, angle_z);
            points.push(p);
        }
    }
    let (_v, _i) = ConvexHull::try_new(&points, None)
        .unwrap()
        .vertices_indices();
}

#[test]
#[ignore]
fn heavy_sphere_test() {
    fn rot_z(point: &Vec<f64>, angle: f64) -> Vec<f64> {
        let e1 = angle.cos() * point[0] - angle.sin() * point[1];
        let e2 = angle.sin() * point[0] + angle.cos() * point[1];
        let e3 = point[2];
        vec![e1, e2, e3]
    }
    fn rot_x(point: &Vec<f64>, angle: f64) -> Vec<f64> {
        let e1 = point[0];
        let e2 = angle.cos() * point[1] - angle.sin() * point[2];
        let e3 = angle.sin() * point[1] + angle.cos() * point[2];
        vec![e1, e2, e3]
    }
    let mut points = Vec::new();
    let dev = 100;
    let unit_y = vec![0.0, 1.0, 0.0];
    for step_x in 0..dev {
        let angle_x = 2.0 * std::f64::consts::PI * (step_x as f64 / dev as f64);
        let p = rot_x(&unit_y, angle_x);
        for step_z in 0..dev {
            let angle_z = 2.0 * std::f64::consts::PI * (step_z as f64 / dev as f64);
            let p = rot_z(&p, angle_z);
            points.push(p);
        }
    }
    let (_v, _i) = ConvexHull::try_new(&points, None)
        .unwrap()
        .vertices_indices();
}

#[test]
fn is_degenerate_test() {
    let points = vec![
        vec![1., 0., 0.],
        vec![0., 0., 0.],
        vec![0., 1., 0.],
        vec![1., 0., 1.],
    ];
    assert!(!is_degenerate(&points));
}
