use ark_ec::{
    models::{short_weierstrass::SWCurveConfig, CurveConfig},
    pairing::{MillerLoopOutput, Pairing, PairingOutput},
    AffineRepr,
};
use ark_ec::CurveGroup;
use ark_ff::{
    fields::{
        fp12_2over3over2::{Fp12, Fp12Config},
        fp2::Fp2Config,
        fp6_3over2::Fp6Config,
        Field, Fp2, PrimeField,
    }, CyclotomicMultSubgroup,
};
use ark_std::{marker::PhantomData};
use ark_ff::{One};
use rayon::prelude::*;


#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// A particular Panther12 group can have G2 being either a multiplicative or a
/// divisive twist.
pub enum TwistType {
    M,
    D,
}

pub trait Panther12Parameters: 'static {
    /// Parameterizes the Panther12 family.
    const T_MINUS_ONE: &'static [u64];

    /// What kind of twist is this?
    const TWIST_TYPE: TwistType;

    // Coefficients required for final expnonentiation
    const C_0:&'static [u64];
    const C_1:&'static [u64];
    const C_2:&'static [u64];
    const C_3:&'static [u64];

    //constant in binary for algorithm 4 of housni gullevic
    const A_0:&'static [u8];
    const A_1:&'static [u8];
    const A_2:&'static [u8];
    const A_3:&'static [u8];

    const LOOP_SIZE: usize;

    const A_NEGATIVE:&'static [bool];

    const MUL_BY_Q_XCOEFF: Fp2<Self::Fp2Config>;
    const MUL_BY_Q_YCOEFF: Fp2<Self::Fp2Config>;

    type Fp: PrimeField + Into<<Self::Fp as PrimeField>::BigInt>;
    type Fp2Config: Fp2Config<Fp = Self::Fp>;
    type Fp6Config: Fp6Config<Fp2Config = Self::Fp2Config>;
    type Fp12Config: Fp12Config<Fp6Config = Self::Fp6Config>;
    type G1Parameters: SWCurveConfig<BaseField = Self::Fp>;
    type G2Parameters: SWCurveConfig<
        BaseField = Fp2<Self::Fp2Config>,
        ScalarField = <Self::G1Parameters as CurveConfig>::ScalarField,
    >;
}

pub mod g1;
pub mod g2;

pub use self::{
    g1::{G1Affine, G1Prepared, G1Projective},
    g2::{G2Affine, G2Prepared, G2Projective, G2HomProjective},
};

#[derive(Derivative)]
#[derivative(Copy, Clone, PartialEq, Eq, Debug, Hash)]
pub struct Panther12<P: Panther12Parameters>(PhantomData<fn() -> P>);

impl<P: Panther12Parameters> Panther12<P> {
    // Evaluate the line function at point p.
    fn ell(f: &mut Fp12<P::Fp12Config>, coeffs: &g2::EllCoeff<P>, p: &G1Affine<P>) {
        let mut c0 = coeffs.0;
        let mut c1 = coeffs.1;
        let mut c2 = coeffs.2;
        let (px, py) = p.xy().unwrap();

        match P::TWIST_TYPE {
            TwistType::M => {
                //multiply f by c_0+(c_2*Py)w+(c_1*Px)vw
                c2.mul_assign_by_fp(py);
                c1.mul_assign_by_fp(px);
                f.mul_by_014(&c0, &c1, &c2);
            },
            TwistType::D => {
                //multiply f by c_0*Py+(c_1*Px)w+(c_2)vw
                c0.mul_assign_by_fp(py);
                c1.mul_assign_by_fp(px);              //                                0       1       2        3             4             5        
                 f.mul_by_034(&c0, &c1, &c2);   //why 034? becuase f is of form (f_0+   f_1v+   f_2v^2+   f_3w+        f_4vw+     f_5v^2w)
                                                      // we multiply it by              c_0*Py+             (c_1*Px)w+      (c_2)vw
            },
        }
    }

   
}

impl<P: Panther12Parameters> Pairing for Panther12<P> {
    type BaseField = <P::G1Parameters as CurveConfig>::BaseField;
    type ScalarField = <P::G1Parameters as CurveConfig>::ScalarField;
    type G1 = G1Projective<P>;
    type G1Affine = G1Affine<P>;
    type G1Prepared = G1Prepared<P>;
    type G2 = G2Projective<P>;
    type G2Affine = G2Affine<P>;
    type G2Prepared = G2Prepared<P>;
    type TargetField = Fp12<P::Fp12Config>;

    fn multi_miller_loop(
        a: impl IntoIterator<Item = impl Into<Self::G1Prepared>>,
        b: impl IntoIterator<Item = impl Into<Self::G2Prepared>>,
    ) -> MillerLoopOutput<Self> {
        use std::ops::Neg;

        // let pairs = a
        //     .into_iter()
        //     .zip_eq(b)
        //     .filter_map(|(p, q)| {
        //         let (p, q) = (p.into(), q.into());
        //         match !p.is_zero() && !q.is_zero() {
        //             true => Some((p, q)),
        //             false => None,
        //         }
        //     })
        //     .collect::<Vec<(Self::G1Prepared,Self::G2Prepared)>>();

        // let f = cfg_chunks_mut!(pairs, 4)
        //     .map(|pairs:Vec<(i64,i64)>| {
        //         let mut f = Self::TargetField::one();
        //         for (p,q) in pairs.iter(){
        //             let pp=q;
        //         }
        //         // for i in BitIteratorBE::without_leading_zeros(P::T_MINUS_ONE).skip(1) {
        //         //     f.square_in_place();
        //         //     for (p, coeffs) in pairs.iter_mut() {
        //         //         Self::ell(&mut f, &coeffs.next().unwrap(), &p.0);
        //         //     }
        //         //     if i {
        //         //         for (p, coeffs) in pairs.iter_mut() {
        //         //             Self::ell(&mut f, &coeffs.next().unwrap(), &p.0);
        //         //         }
        //         //     }
        //         // }
        //         f
        //     })
        //     .product::<Self::TargetField>();

        //algorithm 4 in https://eprint.iacr.org/2021/1359.pdf

        let mut f = Self::TargetField::one();
        for (p,q) in a.into_iter().zip(b.into_iter()).into_iter(){
            let paff = p.into().0;
            let qaff = q.into().0;
            // let Q : [Self::G2Affine;4] =[qaff.clone().neg(),qaff.clone().frobenius(1),qaff.clone().frobenius(2).neg(),qaff.clone().frobenius(3)];

            //instead of \pi_q^i(Q) we compute q^i * Q. These are same becuase in G2, frobenius is same as multiplication by q
            let q1 = mul_by_char::<P>(&qaff);
            let q2 = mul_by_char::<P>(&q1);
            let q3 = mul_by_char::<P>(&q2);

            let mut qvec = [qaff,q1,q2,q3];
            for i in 0..4 {
                if P::A_NEGATIVE[i]{
                    qvec[i]=qvec[i].neg();
                }
            }
            let (t0,t1) = Self::precompute(qvec, &paff);
            let mut j = (P::A_0[0] + 2*P::A_1[0] + 4*P::A_2[0] + 8*P::A_3[0]) as usize;
            let mut f_prime = t1[j-1];
            let s_proj = <Self as Pairing>::G2::from(t0[j-1]);
            let mut s = G2HomProjective::<P>{
                x: s_proj.x,
                y: s_proj.y,
                z: s_proj.z,
            };
            let two_inv = P::Fp::one().double().inverse().unwrap();
            for i in 1..P::LOOP_SIZE{
                f_prime.square_in_place();
                Self::ell(&mut f_prime, &(s.double_in_place(&two_inv)), &paff);
                j = (P::A_0[i] + 2*P::A_1[i] + 4*P::A_2[i] + 8*P::A_3[i]) as usize;

                if j>0 {
                    Self::ell(&mut f_prime,&(s.add_in_place(&t0[j-1])), &paff);
                    if t1[j-1]!=Self::TargetField::one(){
                        f_prime = f_prime*t1[j-1];
                    }
                }
            }
            f = f*f_prime;
        }
        MillerLoopOutput(f)
    }

    fn final_exponentiation(f: MillerLoopOutput<Self>) -> Option<PairingOutput<Self>> {
        // let f=f.0;
        // let q64=P::Fp::characteristic();
        // let r64=P::G2Parameters::get_r();

        // let mut q32=vec!(0 as u32;2*q64.len());
        // let mut r32=vec!(0 as u32;2*r64.len());
        // u64_to_u32(q64, &mut q32[..]);
        // u64_to_u32(r64, &mut r32[..]);

        // let q=BigUint::from_slice(&q32[..]);
        // let r=BigUint::from_slice(&r32[..]);
        // let one=BigUint::from_slice(&[0x1]);
        // let k=12 as u32;

        // let exp=(q.pow(k)-one)/r;
        // let fin=f.pow(&exp.to_u64_digits()[..]);
        // Some(PairingOutput(fin))

        //easy part, f^(q^6-1)
        let mut f=f.0;
        let f_inv=f.inverse().unwrap();
        f = f.frobenius_map(6);
        f=f*f_inv;

        //f^(q^2+1)
        let f_copy=f.clone();
        f = f.frobenius_map(2);
        f=f*f_copy;

        
        // let mut f_3=f.cyclotomic_exp(P::C_3);
        // f_3.frobenius_map(3);
        // let mut f_2=f.cyclotomic_exp(P::C_2);
        // f_2.frobenius_map(2);
        // let mut f_1=f.cyclotomic_exp(P::C_1);
        // f_1.frobenius_map(1);
        // let f_0=f.cyclotomic_exp(P::C_0);

        //hard part f^(c_3 q^3 + c_2 q^2 + c_1 q +c_0)
        let fs=Self::multi_exp_cyclotomic(f, vec![P::C_0,P::C_1,P::C_2,P::C_3]);
        let  f_0=fs[0];
        let mut f_1=fs[1];
        let mut f_2=fs[2];
        let mut f_3=fs[3];
        f_3.frobenius_map_in_place(3);
        f_2.frobenius_map_in_place(2);
        f_1.frobenius_map_in_place(1);
        Some(PairingOutput(f_0*f_1*f_2*f_3))


    }
}

impl<P: Panther12Parameters> Panther12<P>{
    /// compute f*l_{q1,q2}(p)
    fn line(f: Fp12<P::Fp12Config>, q1: G2Affine<P>, q2: G2Affine<P>, p: &G1Affine<P>) -><Self as Pairing>::TargetField{
        let mut f_mut=f.clone();
        let q1_proj = <Self as Pairing>::G2::from(q1);
        let mut q1_proj_hom= G2HomProjective::<P>{
            x: q1_proj.x,
            y: q1_proj.y,
            z: q1_proj.z,
        };
        let ell_coeff = if q1 == q2 {
            let two_inv = P::Fp::one().double().inverse().unwrap();
            q1_proj_hom.double_in_place(&two_inv)
        }
        else{
            q1_proj_hom.add_in_place(&q2)
        };
        Self::ell(&mut f_mut, &ell_coeff, &p);
        f_mut
    }
    ///algorithm 3 in https://eprint.iacr.org/2021/1359.pdf
    fn precompute(q:[G2Affine<P>;4], p: &G1Affine<P>) -> ([G2Affine<P>;15], [Fp12<P::Fp12Config>;15]){
        let pow2=[1 as usize,2 as usize,4 as usize,8 as usize];
        let mut t0=[<Self as Pairing>::G2Affine::zero();15];
        let mut t1=[<Self as Pairing>::TargetField::one();15];
        for i in 0..4{
            t0[pow2[i]-1]=q[i];
        }
        for m in 0..4{
            for n in (m+1)..4{
                let i = pow2[m] + pow2[n];
                t0[i-1]=(t0[pow2[m]-1] + t0[pow2[n]-1]).into_affine();
                let one = <Self as Pairing>::TargetField::one();
                t1[i-1]=Self::line(one,q[m],q[n],&p);
            }
        }

        for m in 0..4{
            for n in (m+1)..4{
                for s in (n+1)..4{
                    let i = pow2[m] + pow2[n] + pow2[s];
                    t0[i-1] = (t0[pow2[m]+pow2[n]-1]+t0[pow2[s]-1]).into_affine();
                    t1[i-1] = Self::line(t1[pow2[m]+pow2[n]-1],(q[m]+q[n]).into_affine(),q[s],&p);
                }
            }
        }
        t0[15-1]=(t0[7-1]+t0[8-1]).into_affine();
        t1[15-1]=Self::line(t1[7-1],(q[0]+q[1]+q[2]).into_affine(),q[3],&p);
        (t0,t1)
    }

    fn multi_exp_cyclotomic(f: Fp12<P::Fp12Config>, v : Vec<&'static [u64]>)->Vec<Fp12<P::Fp12Config>>{
        v.par_iter().map(|e| f.cyclotomic_exp(e)).collect::<Vec<Fp12<P::Fp12Config>>>()
    }
    
    
}

pub fn mul_by_char<P: Panther12Parameters>(r: &G2Affine<P>) -> G2Affine<P> {
    // multiply by field characteristic

    let mut s = *r;
    s.x.frobenius_map_in_place(1);
    s.x *= &P::MUL_BY_Q_XCOEFF;
    s.y.frobenius_map_in_place(1);
    s.y *= &P::MUL_BY_Q_YCOEFF;

    s
}



