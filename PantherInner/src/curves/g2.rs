use crate::fields::{fq::*, fr::*, fq2::*};
use ark_ec::{
    models::CurveConfig,
    short_weierstrass::{self}
};
use ark_ff::{ MontFp, Zero};
use crate::panther12;

pub type G2Affine = panther12::G2Affine<crate::curves::pairing::Parameters>;
pub type G2Projective = panther12::G2Projective<crate::curves::pairing::Parameters>;

#[derive(Clone, Default, PartialEq, Eq)]
pub struct Parameters;

impl CurveConfig for Parameters {
    type BaseField = Fq2;
    type ScalarField = Fr;

    
    #[rustfmt::skip]
    const COFACTOR: &'static [u64] = &[
        0xc48c51435e50d794,
        0xc9c2d8269401c38b,
        0xbc31fc96ed25063d,
        0x648ac7157eeaa4ef,
        0x47f1ca60d809e6bc,
        0xd655f63bd9839262,
        0x38b41c71cc26d95e,
        0x327ab8666c75cd54,
        0xdfa1fb797e3c42e0,
        0xc735c81d2e71516a,
        0x3fa21bb14969c708,
        0x5985f2aec98703b6,
        0x2b667470f4d28ff5,
        0x73d0615cf5dceef4,
        0x9d3a1acd512,
    ];
    /// COFACTOR_INV = COFACTOR^{-1} mod r
    /// 50643695498298099292022683002431747179445631306923841235773964100107627076818
    #[rustfmt::skip]
    const COFACTOR_INV: Fr = MontFp!(
        "50643695498298099292022683002431747179445631306923841235773964100107627076818"
    );
}

impl short_weierstrass::SWCurveConfig for Parameters {
    /// COEFF_A = [0, 0]
    const COEFF_A: Fq2 = Fq2::new(MontFp!("0"), MontFp!("0"));
    //we are using D twist
    /// COEFF_B = 14/(7u+49)
    const COEFF_B: Fq2 = Fq2::new(
                                    MontFp!("383243406726075367023340369723164497849297580969953147685129555230975819662660710383675930554440305511012701911036594658199726468011940153311225682737268035151061272723571376914433"),
                                    MontFp!("355868877674212840807387486171509890860062039472099351421906015571620403972470659641984792657694569403083223203105409325471174577439658713788995276827463175497414038957601992849116")
                                );

    /// AFFINE_GENERATOR_COEFFS = (G2_GENERATOR_X, G2_GENERATOR_Y)
    const GENERATOR: G2Affine = G2Affine::new_unchecked(G2_GENERATOR_X, G2_GENERATOR_Y);

    #[inline(always)]
    fn mul_by_a(_: Self::BaseField) -> Self::BaseField {
        Self::BaseField::zero()
    }

    // fn is_in_correct_subgroup_assuming_on_curve(point: &G2Affine) -> bool {
    //     // Algorithm from Section 4 of https://eprint.iacr.org/2021/1130.
    //     //
    //     // Checks that [p]P = [T_MINUS_ONE]P
    //     let t1=crate::panther_inner::Parameters::T_MINUS_ONE;
    //     let x_times_point =
    //         point.mul_bigint(BigInt::new([t1[0],t1[1],t1[2],t1[3],t1[4]]));

    //     //let p_times_point = p_power_endomorphism(point);
    //     let p_times_point = point.mul_bigint(FqConfig::MODULUS);

    //     x_times_point.eq(&p_times_point)
    // }

    
}

pub const G2_GENERATOR_X: Fq2 = Fq2::new(G2_GENERATOR_X_C0, G2_GENERATOR_X_C1);
pub const G2_GENERATOR_Y: Fq2 = Fq2::new(G2_GENERATOR_Y_C0, G2_GENERATOR_Y_C1);

/// G2_GENERATOR_X_C0 =
/// 352701069587466618187139116011060144890029952792775240219908644239793785735715026873347600343865175952761926303160
#[rustfmt::skip]
pub const G2_GENERATOR_X_C0: Fq = MontFp!("132247468448863276600327856158541757108730546311256563312960653795693648960560183198026253798999106583125159920675806114964109122303030229678263356594223243024175648188692044868140");

/// G2_GENERATOR_X_C1 =
/// 3059144344244213709971259814753781636986470325476647558659373206291635324768958432433509563104347017837885763365758
#[rustfmt::skip]
pub const G2_GENERATOR_X_C1: Fq = MontFp!("137757259767192542734174498153039462622956842398003407509988598850157274712372253550878341178160392516522592233382999402508197531853239295683360903545197489829784833018300080781707");

/// G2_GENERATOR_Y_C0 =
/// 1985150602287291935568054521177171638300868978215655730859378665066344726373823718423869104263333984641494340347905
#[rustfmt::skip]
pub const G2_GENERATOR_Y_C0: Fq = MontFp!("398378319740754339326337869079960873246804305549431111500421860262585598577930300582421299597240289501650483784725698636583897231884929044709341582778823859388177814822672039221246");

/// G2_GENERATOR_Y_C1 =
/// 927553665492332455747201965776037880757740193453592970025027978793976877002675564980949289727957565575433344219582
#[rustfmt::skip]
pub const G2_GENERATOR_Y_C1: Fq = MontFp!("289144805579014228712235364663577703036634319484951854067827214830103555510371544251320903950673842086416059171530835034445879364767726631330236201856913107981838247409575758434545");



#[test]
fn test_curve(){
    use ark_ff::PrimeField;
    use ark_ec::AffineRepr;
    use ark_ff::UniformRand;
    use num_bigint::BigUint;
    use ark_ff::Field;
    use super::util::u64_to_u32;
    use crate::panther12::mul_by_char;
    use super::pairing::Parameters as P;
    use super::g2::Parameters as PG2;
    use ark_ec::CurveGroup;
    let scalar= <PG2 as CurveConfig>::ScalarField::MODULUS;
    let q = <PG2 as CurveConfig>::BaseField::characteristic();
    let rng = &mut ark_std::test_rng();
    let point = G2Affine::rand(rng);

    //test mul_by_char
    assert_eq!(point.mul_bigint(q).into_affine(),mul_by_char::<P>(&point));

    //check r*generator is infinity
    assert!(point.mul_bigint(scalar).is_zero());

    let scalar=BigUint::from(scalar);
    let cofactor64=<PG2 as CurveConfig>::COFACTOR;
    let mut cofactor32=vec!(0 as u32;2*cofactor64.len());
    u64_to_u32(cofactor64, &mut cofactor32[..]);

    let cofactor=BigUint::from_slice(&cofactor32[..]);

    //check r*cofactor=cardinality
    let scalar_slice:&[u32]=&[4294967277, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 2147483647];
    let cofactor_slice:&[u32]=&[1582356372, 3297530179, 2483143563, 3384989734, 3978626621, 3157392534, 2129306863, 1686816533, 3624527548, 1207028320, 3649278562, 3595957819, 3425098078, 951327857, 1819659604, 846903398, 2117878496, 3751934841, 779178346, 3342190621, 1231669000, 1067588529, 3381068726, 1501950638, 4107440117, 728134768, 4124897012, 1943036252, 2712458514, 2515];
    let cardinality_slice:&[u32]=&[4, 1771436032, 64912544, 109704483, 1715505514, 139083980, 2492842549, 2310224231, 2792115150, 2337547121, 624755834, 4236256674, 1336959252, 2830819801, 2998342967, 1932112880, 227278010, 2331196204, 2052636226, 2720866547, 1933158152, 1666318352, 1094033443, 1951160634, 2474471994, 926308698, 1473440918, 3408112168, 2766213925, 533746467, 1690534363, 2898458967, 2053720058, 364067384, 2062448506, 971518126, 3503712905, 1257];
    assert_eq!(scalar,BigUint::from_slice(scalar_slice));
    assert_eq!(cofactor,BigUint::from_slice(cofactor_slice));
    assert_eq!(cofactor*scalar,BigUint::from_slice(cardinality_slice)); 
}



