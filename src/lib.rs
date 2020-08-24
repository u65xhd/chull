//! Convex hull approximation for rust based on [Quick hull](http://citeseerx.ist.psu.edu/viewdoc/summary;jsessionid=C57E2269B0D64504B97E8469F6A1315D?doi=10.1.1.117.405).
//! Available in 3-D or 2-D for now.
//! ## Examples
//! ```
//! use chull::ConvexHull;
//!
//! let p1 = vec![1.0, 1.0, 1.0];
//! let p2 = vec![1.0, 1.0, -1.0];
//! let p3 = vec![1.0, -1.0, 1.0];
//! let p4 = vec![1.0, -1.0, -1.0];
//! let p5 = vec![-1.0, 1.0, 1.0];
//! let p6 = vec![-1.0, 1.0, -1.0];
//! let p7 = vec![-1.0, -1.0, 1.0];
//! let p8 = vec![-1.0, -1.0, -1.0];
//! let p9 = vec![0.0, 0.0, 0.0];
//! // threshold is used internally to determine whether a point is above or below the surface.
//! let threshold = 0.001;
//! let points = vec![p1, p2, p3, p4, p5, p6, p7, p8, p9];
//! let cube = ConvexHull::try_new(&points, threshold, None).unwrap();
//! assert_eq!(cube.volume(), 8.0);
//! let (_v,i) = cube.vertices_indices();
//! assert_eq!(i.len(), 6 * 2 * 3);
//! ```
//! ## If the results are inaccurate
//! If the calculation results are inaccurate due to rounding errors, an error may occur. In such cases, the use of integer types such as [```i128```](https://doc.rust-lang.org/std/primitive.i128.html) or [```BigDecimal```](https://docs.rs/bigdecimal/) may improve the result.
//! ```
//! use chull::ConvexHull;
//!
//! let p1 = vec![1.0, 0.0, 0.0];
//! let p2 = vec![0.0, 0.001, 0.0];
//! let p3 = vec![0.0, 0.0, 0.00001];
//! let p4 = vec![-1.0, 0.0, 0.0];
//! let p5 = vec![0.0, -0.001, 0.0];
//! let p6 = vec![0.0, 0.0, -0.00001];
//! let points_float = vec![p1, p2, p3, p4, p5, p6];
//! let mut points_i128 = Vec::new();
//! for point_float in &points_float{
//!     points_i128.push(point_float.iter().map(|x| (x*1_000_000.0) as i128).collect::<Vec<_>>());
//! }
//! let octahedron = ConvexHull::try_new(&points_i128, 10, None).unwrap();
//! // The following are likely to be errors
//! //let octahedron = ConvexHull::try_new(&points_float, std::f32::EPSILON, None).unwrap();
//! ```

pub mod convex;
pub mod util;
pub use convex::*;
