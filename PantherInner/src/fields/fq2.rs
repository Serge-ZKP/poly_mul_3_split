use crate::fields::fq::*;
use ark_ff::{fields::*, MontFp};

pub type Fq2 = Fp2<Fq2Config>;

pub struct Fq2Config;

impl Fp2Config for Fq2Config {
    type Fp = Fq;

    /// NONRESIDUE = 7
    #[rustfmt::skip]
    const NONRESIDUE: Fq = MontFp!("7");   //quotienting by u^2-7

    /// Coefficients for the Frobenius automorphism.
    #[rustfmt::skip]
    const FROBENIUS_COEFF_FP2_C1: &'static [Fq] = &[
        // Fq(7)**(((q^0) - 1) / 2)
        MontFp!("1"),
        // Fq(7)**(((q^1) - 1) / 2)
        MontFp!("-1"),
    ];
    
}

pub const FQ2_ZERO: Fq2 = Fq2::new(FQ_ZERO, FQ_ZERO);
pub const FQ2_ONE: Fq2 = Fq2::new(FQ_ONE, FQ_ZERO);

#[test]
fn test_fq2(){
    //check u^2=7
    let u=Fq2::new(MontFp!("0"),MontFp!("1"));
    let seven=Fq2::new(MontFp!("7"),MontFp!("0"));
    assert_eq!(u.square(),seven);
}

#[test]
pub fn test_frobenius() {
    use ark_ff::Field;
    use ark_std::UniformRand;
    let mut rng = ark_std::test_rng();
    let characteristic = Fq2::characteristic();
    let max_power = (Fq2::extension_degree() + 1) as usize;

    for _ in 0..10 {
        let a = Fq2::rand(&mut rng);

        let mut a_0 = a;
        a_0 = a_0.frobenius_map(0);
        assert_eq!(a, a_0);

        let mut a_q = a.pow(&characteristic);
        for power in 1..max_power {
            let mut a_qi = a;
            a_qi = a_qi.frobenius_map(power);
            assert_eq!(a_qi, a_q, "failed on power {}", power);

            a_q = a_q.pow(&characteristic);
        }
    }
}
