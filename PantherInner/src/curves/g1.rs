use ark_ec::{
    models::CurveConfig,
    short_weierstrass::{self, *},
};
use ark_ff::{MontFp, Zero};

use crate::fields::{fq::*, fr::*};

pub type G1Affine = Affine<Parameters>;
pub type G1Projective = Projective<Parameters>;

#[derive(Clone, Default, PartialEq, Eq)]
pub struct Parameters;

impl CurveConfig for Parameters {
    type BaseField = Fq;
    type ScalarField = Fr;
    const COFACTOR: &'static [u64] = &[0x84651050d79435e5, 0x2670092267680b9b,0x48a957910ddd0e5c,0xa6549f829e7cd99a,0x96cf524ce975dea,0x46ee70];
    const COFACTOR_INV: Fr = MontFp!("55054288129687180500829749807481096807670042311671699761348917096905673330423");
}

impl short_weierstrass::SWCurveConfig for Parameters {
    //using curve y^2=x^3+14

    /// COEFF_A = 0
    const COEFF_A: Fq = MontFp!("0");

    /// COEFF_B = 14
    const COEFF_B: Fq = MontFp!("14");

    /// AFFINE_GENERATOR_COEFFS = (G1_GENERATOR_X, G1_GENERATOR_Y)
    const GENERATOR: G1Affine = G1Affine::new_unchecked(G1_GENERATOR_X, G1_GENERATOR_Y);

    #[inline(always)]
    fn mul_by_a(_: Self::BaseField) -> Self::BaseField {
        Self::BaseField::zero()
    }

    
}


/// G1_GENERATOR_X =
/// 266366668494810926016963514089155757623307897246933888951043686670988276065532742878536193527711862671017608170386674119249463021638514914545726765972511315843515110224210578390813
#[rustfmt::skip]
pub const G1_GENERATOR_X: Fq = MontFp!("266366668494810926016963514089155757623307897246933888951043686670988276065532742878536193527711862671017608170386674119249463021638514914545726765972511315843515110224210578390813");

/// G1_GENERATOR_Y =
/// 532483566563071017598378510349051009611057946307036137582731679855439436884875246014054184174950988090164326837880753973682919218152879675656914147764919331721967507643061430234293
#[rustfmt::skip]
pub const G1_GENERATOR_Y: Fq = MontFp!("-532483566563071017598378510349051009611057946307036137582731679855439436884875246014054184174950988090164326837880753973682919218152879675656914147764919331721967507643061430234293");

#[cfg(test)]
mod test {

    use super::*;
    use ark_ec::CurveGroup;
    use ark_std::UniformRand;

    #[test]
    fn batch_normalization() {
        let mut rng = ark_std::test_rng();

        let mut g_s = [G1Projective::zero(); 100];
        for i in 0..100 {
            g_s[i] = G1Projective::rand(&mut rng);
        }

        let mut g_s_affine_naive = [G1Affine::identity(); 100];
        for (i, g) in g_s.iter().enumerate() {
            g_s_affine_naive[i] = g.into_affine();
        }

        let g_s_affine_fast = G1Projective::normalize_batch(&g_s);
        assert_eq!(g_s_affine_naive.as_ref(), g_s_affine_fast.as_slice());
    }

    #[test]
    fn test_curve(){
        use ark_ff::{PrimeField};
        use ark_ec::AffineRepr;
        use num_bigint::BigUint;
        use self::Parameters as PG1;
        use super::super::util::*;
        let scalar= <PG1 as CurveConfig>::ScalarField::MODULUS;
        let point = <PG1 as short_weierstrass::SWCurveConfig>::GENERATOR; 

        //check r*generator is infinity
        assert!(point.mul_bigint(scalar).is_zero());

        let scalar=BigUint::from(scalar);
        let cofactor64=<PG1 as CurveConfig>::COFACTOR;
        let mut cofactor32=vec!(0 as u32;2*cofactor64.len());
        u64_to_u32(cofactor64, &mut cofactor32[..]);

        let cofactor=BigUint::from_slice(&cofactor32[..]);

        //check r*cofactor=cardinality
        let scalar_slice:&[u32]=&[4294967277, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 2147483647];
        let cofactor_slice:&[u32]=&[3616814565, 2221215824, 1734871963, 644876578, 232590940, 1219057553, 2658982298, 2790563714, 3466026474, 158135588, 4648560];
        let cardinality_slice:&[u32]=&[1, 746572288, 1397171061, 632246898, 4170706729, 2607710267, 1018943884, 666380622, 378413704, 253515373, 779113341, 322438289, 2263779118, 609528776, 1329491149, 1395281857, 1733013237, 79067794, 2324280];
        assert_eq!(scalar,BigUint::from_slice(scalar_slice));
        assert_eq!(cofactor,BigUint::from_slice(cofactor_slice));
        assert_eq!(cofactor*scalar,BigUint::from_slice(cardinality_slice));

        
    }

    #[test]
    fn test_cofactor_ops() {
        use ark_ec::AffineRepr;
        let rng = &mut ark_std::test_rng();
        for _ in 0..10{
            let a = G1Affine::rand(rng);
            assert_eq!(a.mul_by_cofactor_to_group(), a.mul_bigint(<super::super::g1::Parameters as CurveConfig>::COFACTOR));
            assert_eq!(a.mul_by_cofactor(), a.mul_bigint(<super::super::g1::Parameters as CurveConfig>::COFACTOR));
            assert_eq!(a.mul_by_cofactor().mul_by_cofactor_inv(), a);
            assert_eq!(a.mul_by_cofactor_inv().mul_by_cofactor(), a);
            assert_eq!(a.mul_by_cofactor_inv(), a * <super::super::g1::Parameters as CurveConfig>::COFACTOR_INV);

            assert!(a.clear_cofactor().is_in_correct_subgroup_assuming_on_curve());
        }
    }
}
