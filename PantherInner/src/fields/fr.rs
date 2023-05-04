use ark_ff::fields::{Fp256, MontBackend, MontConfig};

#[derive(MontConfig)]
#[modulus = "57896044618658097711785492504343953926634992332820282019728792003956564819949"]
#[generator = "2"]
pub struct FrConfig;
pub type Fr = Fp256<MontBackend<FrConfig, 4>>;

#[test]
fn test_inv() {
    assert_eq!(FrConfig::INV, 0x86bc_a1af_286b_ca1b);  // -modulus^{-1} mod 2^64
}

#[test]
fn test_modulus() {
    assert_eq!(
        FrConfig::MODULUS.0,
        [
            0xffff_ffff_ffff_ffed,
            0xffff_ffff_ffff_ffff,
            0xffff_ffff_ffff_ffff,
            0x7fff_ffff_ffff_ffff,
        ]
    );
}

#[test]
pub fn test_frobenius() {
    use ark_ff::Field;
    use ark_std::UniformRand;
    let mut rng = ark_std::test_rng();
    let characteristic = Fr::characteristic();
    let max_power = (Fr::extension_degree() + 1) as usize;

    for _ in 0..10 {
        let a = Fr::rand(&mut rng);
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
