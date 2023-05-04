use ark_ec::{short_weierstrass::{Affine, Projective}, AffineRepr, CurveGroup};
use ark_ff::fields::{Field, Fp2};
use ark_serialize::*;
use ark_std::vec::Vec;
use ark_ec::short_weierstrass::SWCurveConfig;

use crate::panther12::{Panther12Parameters, TwistType};

pub type G2Affine<P> = Affine<<P as Panther12Parameters>::G2Parameters>;
pub type G2Projective<P> = Projective<<P as Panther12Parameters>::G2Parameters>;

#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(
    Clone(bound = "P: Panther12Parameters"),
    Debug(bound = "P: Panther12Parameters"),
    PartialEq(bound = "P: Panther12Parameters"),
    Eq(bound = "P: Panther12Parameters")
)]
pub struct G2Prepared<P: Panther12Parameters>(pub G2Affine<P>);

pub(crate) type EllCoeff<P> = (
    Fp2<<P as Panther12Parameters>::Fp2Config>,
    Fp2<<P as Panther12Parameters>::Fp2Config>,
    Fp2<<P as Panther12Parameters>::Fp2Config>,
);

#[derive(Derivative)]
#[derivative(
    Clone(bound = "P: Panther12Parameters"),
    Copy(bound = "P: Panther12Parameters"),
    Debug(bound = "P: Panther12Parameters")
)]
pub struct G2HomProjective<P: Panther12Parameters> {
    pub x: Fp2<P::Fp2Config>,
    pub y: Fp2<P::Fp2Config>,
    pub z: Fp2<P::Fp2Config>,
}

impl<P: Panther12Parameters> Default for G2Prepared<P> {
    fn default() -> Self {
        Self::from(G2Affine::<P>::generator())
    }
}

impl<P: Panther12Parameters> From<G2Affine<P>> for G2Prepared<P> {
    fn from(q: G2Affine<P>) -> Self {
        G2Prepared(q)
    }
}

impl<P: Panther12Parameters> G2Prepared<P> {
    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}


impl<P: Panther12Parameters> From<G2Projective<P>> for G2Prepared<P> {
    fn from(q: G2Projective<P>) -> Self {
        q.into_affine().into()
    }
}

impl<'a, P: Panther12Parameters> From<&'a G2Affine<P>> for G2Prepared<P> {
    fn from(other: &'a G2Affine<P>) -> Self {
        (*other).into()
    }
}

impl<'a, P: Panther12Parameters> From<&'a G2Projective<P>> for G2Prepared<P> {
    fn from(q: &'a G2Projective<P>) -> Self {
        q.into_affine().into()
    }
}

impl<P: Panther12Parameters> G2HomProjective<P> {
    pub fn double_in_place(&mut self, two_inv: &P::Fp) -> EllCoeff<P> {

        //double in place and return the line function

        // Formula for line function when working with
        // homogeneous projective coordinates.

        let mut a = self.x * &self.y;
        a.mul_assign_by_fp(two_inv);
        let b = self.y.square();
        let c = self.z.square();
        let e = P::G2Parameters::COEFF_B * &(c.double() + &c);
        let f = e.double() + &e;
        let mut g = b + &f;
        g.mul_assign_by_fp(two_inv);
        let h = (self.y + &self.z).square() - &(b + &c);
        let i = e - &b;
        let j = self.x.square();
        let e_square = e.square();

        self.x = a * &(b - &f);
        self.y = g.square() - &(e_square.double() + &e_square);
        self.z = b * &h;
        match P::TWIST_TYPE {
            TwistType::M => (i, j.double() + &j, -h),
            TwistType::D => (-h, j.double() + &j, i),
        }
    }

    pub fn add_in_place(&mut self, q: &G2Affine<P>) -> EllCoeff<P> {
        let (&qx, &qy) = q.xy().unwrap();
        // add the point and return the line function
        
        // Formula for line function when working with
        // homogeneous projective coordinates.
        let theta = self.y - &(qy * &self.z);
        let lambda = self.x - &(qx * &self.z);
        let c = theta.square();
        let d = lambda.square();
        let e = lambda * &d;
        let f = self.z * &c;
        let g = self.x * &d;
        let h = e + &f - &g.double();
        self.x = lambda * &h;
        self.y = theta * &(g - &h) - &(e * &self.y);
        self.z *= &e;
        let j = theta * &qx - &(lambda * &qy);

        match P::TWIST_TYPE {
            TwistType::M => (j, -theta, lambda),
            TwistType::D => (lambda, -theta, j),
        }
    }
}

