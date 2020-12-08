# chull
![Build Status](https://github.com/daidanretsu/chull/workflows/Rust/badge.svg) [![crates.io](https://img.shields.io/crates/v/chull.svg)](https://crates.io/crates/chull)

[Documentation](https://docs.rs/chull/)
![](img/usagi.png)
Convex hull approximation based on [Quick hull](http://citeseerx.ist.psu.edu/viewdoc/summary;jsessionid=C57E2269B0D64504B97E8469F6A1315D?doi=10.1.1.117.405).
Available in 3-D or 2-D for now.
## Examples
```rust
use chull::{ConvexHull,ConvexHullWrapper};
let p1 = vec![1.0, 1.0, 1.0];
let p2 = vec![1.0, 1.0, -1.0];
let p3 = vec![1.0, -1.0, 1.0];
let p4 = vec![1.0, -1.0, -1.0];
let p5 = vec![-1.0, 1.0, 1.0];
let p6 = vec![-1.0, 1.0, -1.0];
let p7 = vec![-1.0, -1.0, 1.0];
let p8 = vec![-1.0, -1.0, -1.0];
let p9 = vec![0.0, 0.0, 0.0];

let points = vec![p1, p2, p3, p4, p5, p6, p7, p8, p9];

// slow but robust
let cube1 = ConvexHullWrapper::try_new(&points, None).unwrap();
// fast but non-robust
let cube2 = ConvexHull::try_new(&points, std::f64::EPSILON*200.0, None).unwrap();

assert_eq!(cube1.volume(), 8.0);
let (v,i) = cube1.vertices_indices();
assert_eq!(v.len(), 8);
assert_eq!(i.len(), 6 * 2 * 3);
```

License: MIT/Apache-2.0
