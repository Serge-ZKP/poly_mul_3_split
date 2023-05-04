use crate::panther12::Panther12Parameters;
use ark_ec::{short_weierstrass::{Affine, Projective}, CurveGroup, AffineRepr};
use ark_serialize::*;
use ark_std::vec::Vec;

pub type G1Affine<P> = Affine<<P as Panther12Parameters>::G1Parameters>;
pub type G1Projective<P> = Projective<<P as Panther12Parameters>::G1Parameters>;

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(
    Clone(bound = "P: Panther12Parameters"),
    Debug(bound = "P: Panther12Parameters"),
    PartialEq(bound = "P: Panther12Parameters"),
    Eq(bound = "P: Panther12Parameters")
)]
pub struct G1Prepared<P: Panther12Parameters>(pub G1Affine<P>);

impl<P: Panther12Parameters> From<G1Affine<P>> for G1Prepared<P> {
    fn from(other: G1Affine<P>) -> Self {
        G1Prepared(other)
    }
}

impl<P: Panther12Parameters> From<G1Projective<P>> for G1Prepared<P> {
    fn from(q: G1Projective<P>) -> Self {
        q.into_affine().into()
    }
}

impl<'a, P: Panther12Parameters> From<&'a G1Affine<P>> for G1Prepared<P> {
    fn from(other: &'a G1Affine<P>) -> Self {
        G1Prepared(*other)
    }
}

impl<'a, P: Panther12Parameters> From<&'a G1Projective<P>> for G1Prepared<P> {
    fn from(q: &'a G1Projective<P>) -> Self {
        q.into_affine().into()
    }
}

impl<P: Panther12Parameters> G1Prepared<P> {
    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl<P: Panther12Parameters> Default for G1Prepared<P> {
    fn default() -> Self {
        G1Prepared(G1Affine::<P>::generator())
    }
}
