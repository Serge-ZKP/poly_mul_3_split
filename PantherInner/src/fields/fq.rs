use ark_ff::fields::{Fp640, MontBackend};

#[derive(ark_ff::MontConfig)]
#[modulus = "574865110089113050535010554584746746773946371454929721527694332846463729493991065575513895831660458266519052866554891987299589702017910229966838524105902052726591909085357065371649"]
#[generator = "14"]
pub struct FqConfig;
pub type Fq = Fp640<MontBackend<FqConfig, 10>>;  //10 is the number of 64 bit integers required to represent an element of Fq


pub const FQ_ONE: Fq = ark_ff::MontFp!("1");
pub const FQ_ZERO: Fq = ark_ff::MontFp!("0");


#[test]
pub fn test_frobenius() {
    use ark_ff::Field;
    use ark_std::UniformRand;
    let mut rng = ark_std::test_rng();
    let characteristic = Fq::characteristic();
    let max_power = (Fq::extension_degree() + 1) as usize;

    for _ in 0..10 {
        let a = Fq::rand(&mut rng);

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
