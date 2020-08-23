use std::cmp::Ordering;
use std::ops::{Add,Sub,Mul,Div,Rem};
use bigdecimal::{BigDecimal,ParseBigDecimalError};
use num_bigint::BigInt;
use num_traits::{Zero,One,Num};
use super::traits::*;

//impl Sqrt for BigDecimal {
//    fn sqrt(&self) -> Option<BigDecimal> {
//        self.sqrt()
//    }
//}
//
//impl ApproxSign for BigDecimal {
//    fn approx_sign(&self) -> Sign {
//        let threshold = BigDecimal::new(BigInt::new(num_bigint::Sign::Plus, vec![1]), 20);
//        if self.abs() <= threshold {
//            Sign::NoSign
//        } else if self > &threshold {
//            Sign::Plus
//        } else {
//            Sign::Minus
//        }
//    }
//}

#[derive(Clone)]
pub (crate) struct BDWrapper(pub (crate) BigDecimal);

impl Sqrt for BDWrapper {
    fn sqrt(&self) -> Option<BDWrapper> {
        self.0.sqrt().map(|s| Self(s))
    }
}


impl PartialEq for BDWrapper {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}


impl PartialOrd for BDWrapper {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl Add for BDWrapper {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self (self.0+other.0)
    }
}

impl Sub for BDWrapper {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self (self.0-other.0)

    }
}

impl Mul for BDWrapper {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Self (self.0*rhs.0)

    }
}

impl Div for BDWrapper {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        Self (self.0/rhs.0)

    }
}

impl Rem for BDWrapper {
    type Output = Self;

    fn rem(self, modulus: Self) -> Self {
        Self (self.0%modulus.0)

    }
}

impl Zero for BDWrapper{
    fn zero() -> Self{
        Self(BigDecimal::zero())
    }
    fn is_zero(&self) -> bool{
        self.0.is_zero()
    }
}

impl One for BDWrapper{
    fn one() -> Self{
        Self(BigDecimal::one())
    }
}

impl Num for BDWrapper{
    type FromStrRadixErr = ParseBigDecimalError;
    fn from_str_radix(str: &str, radix: u32) -> Result<Self, Self::FromStrRadixErr>
    {
        BigDecimal::from_str_radix(str, radix).map(|bd| Self(bd))
    }
}

impl ApproxSign for BDWrapper {
    fn approx_sign(&self) -> Sign {
        let threshold = BigDecimal::new(BigInt::new(num_bigint::Sign::Plus, vec![1]), 20);
        if self.0.abs() <= threshold {
            Sign::NoSign
        } else if self.0 > threshold {
            Sign::Plus
        } else {
            Sign::Minus
        }
    }
}