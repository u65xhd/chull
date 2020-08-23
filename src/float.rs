use num_traits::{One,Float};
use super::traits::*;

impl<T: Float> Sqrt for T 
{
    fn sqrt(&self) -> Option<Self> {
        let sqrt = self.clone().sqrt();
        if sqrt.is_nan() {
            None
        } else {
            Some(sqrt)
        }
    }
}

impl<T: Float> ApproxSign for T {
    fn approx_sign(&self) -> Sign {
        let ten = T::one()+T::one()+T::one()+T::one()+T::one()+T::one()+T::one()+T::one()+T::one()+T::one();
        let threshold = T::epsilon();// * ten*ten+(T::one()+T::one());
        if self.abs() <= threshold {
            Sign::NoSign
        } else if self > &threshold {
            Sign::Plus
        } else {
            Sign::Minus
        }
    }
}