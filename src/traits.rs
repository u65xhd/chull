use std::marker::Sized;

pub trait Sqrt {
    fn sqrt(&self) -> Option<Self>
    where
        Self: Sized;
}

#[derive(Debug, Clone, Copy, PartialOrd, PartialEq, Eq, Hash)]
pub enum Sign {
    Minus,
    NoSign,
    Plus,
}

pub trait ApproxSign {
    fn approx_sign(&self) -> Sign
    where
        Self: Sized;
}