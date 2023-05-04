
    #[test]
    fn test_bilinearity() {
        use ark_ec::{pairing::*};
        use ark_std::{test_rng, UniformRand};
        use ark_ff::Field;
        use ark_ff::Zero;
        use ark_ec::Group;
        let mut rng = test_rng();
        let a: <crate::curves::Panther12_598 as Pairing>::G1 = UniformRand::rand(&mut rng);
        let b: <crate::curves::Panther12_598 as Pairing>::G2 = UniformRand::rand(&mut rng);
        let s: <crate::curves::Panther12_598 as Pairing>::ScalarField = UniformRand::rand(&mut rng);

        let sa = a * s;
        let sb = b * s;

        let ans1 = <crate::curves::Panther12_598>::pairing(sa, b);
        let ans2 = <crate::curves::Panther12_598>::pairing(a, sb);
        let ans3 = <crate::curves::Panther12_598>::pairing(a, b) * s;


        assert_eq!(ans1.0, ans2.0);
        assert_eq!(ans2, ans3);

        assert_ne!(ans1, PairingOutput::zero());
        assert_ne!(ans2, PairingOutput::zero());
        assert_ne!(ans3, PairingOutput::zero());
        let group_order = <<crate::curves::Panther12_598 as Pairing>::ScalarField>::characteristic();

        assert_eq!(ans1.mul_bigint(group_order), PairingOutput::zero());
        assert_eq!(ans2.mul_bigint(group_order), PairingOutput::zero());
        assert_eq!(ans3.mul_bigint(group_order), PairingOutput::zero());
    }

    #[test]
    fn time_pairing() {
        use ark_ec::short_weierstrass::SWCurveConfig;
        use ark_ec::AffineRepr;
        use ark_ec::{pairing::*};
        use crate::curves::g1::Parameters as PG1;
        use crate::curves::g2::Parameters as PG2;
        use std::time::Instant;
        let now = Instant::now();
        let a: <crate::curves::Panther12_598 as Pairing>::G1 = PG1::GENERATOR.into_group();
        let b: <crate::curves::Panther12_598 as Pairing>::G2 = PG2::GENERATOR.into_group();
        let _c = <crate::curves::Panther12_598>::pairing(a, b);
        let elapsed = now.elapsed();
        println!("Elapsed: {:.2?}", elapsed);
    }

    #[test]
    fn test_msm() {
        use crate::curves::g1::{G1Projective as G, G1Affine as GAffine};
        use ark_ff::UniformRand;
        use std::time::Instant;
        use crate::fields::fr::Fr as ScalarField;
        use ark_ec::VariableBaseMSM;
        let mut rng = ark_std::test_rng();
        // Let's sample uniformly random group elements:
        const DEG: usize = 1000000;
        let mut g = vec![];
        let mut s = vec![];
        for _ in 0..DEG {
            g.push(GAffine::rand(&mut rng));
            s.push(ScalarField::rand(&mut rng));
        }

        let now = Instant::now();
        let _r = G::msm(g.as_slice(), s.as_slice()).unwrap();
        let elapsed = now.elapsed();
        println!("Elapsed: {:.2?}", elapsed);
    }

