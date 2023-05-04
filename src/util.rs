use std::cmp::min;

use ark_ec::{pairing::Pairing, AffineRepr};
use ark_ff::{PrimeField};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};
use ark_std::rand::rngs::StdRng;
use rand_core::SeedableRng;
use ark_ec::CurveGroup;
use rayon::{iter::IntoParallelRefIterator, prelude::{IndexedParallelIterator, ParallelIterator}};

use crate::wrapper_func::mul_ark_poly;

/// Prepare the CRS and verifer key. Used only for testing purposes
#[macro_export]
macro_rules! prepare_kzg10 {
    // `()` indicates that the macro takes no argument.
    ($N: expr, $powers: ident, $vk: ident) => {
        use crate::commitment::{ExtendedCommitment, ExtendedKZG10};
        use panther_inner::curves::Panther12_598;
        use ark_std::test_rng;
        type KZG= ExtendedKZG10<Panther12_598>;

        let rng = &mut test_rng();
        let pp = KZG::extended_setup(4*$N, rng);
        let ($powers, $vk) = KZG::trim(&pp, $N *2, $N);
}
}

/// used for logging
pub fn init_log() {
    let _ = env_logger::builder().is_test(true).try_init();
}

/// create random polynonmial of degree `deg`
pub fn rand_poly<F: PrimeField>(deg: usize) -> DensePolynomial<F>{
    let rng = &mut StdRng::from_entropy();
    let mut coeffs = vec![];
    for _ in 0..deg+1{
        coeffs.push(F::rand(rng));
    }
    DensePolynomial::from_coefficients_vec(coeffs)
}

#[allow(non_snake_case)]
/// return `x^n * f(1/x)`
pub fn x_power_n_times_f_of_x_inverse<F: PrimeField>(f: &DensePolynomial<F>, N: usize) -> DensePolynomial<F>{
    let mut coeffs = f.coeffs().to_vec();
    while coeffs.len() <= N{
        coeffs.push(F::zero());
    }
    DensePolynomial::<F>::from_coefficients_vec(coeffs.iter().copied().rev().collect())
}
#[allow(non_snake_case)]
/// return `f(gamma * x)`
pub fn f_gamma_x<F: PrimeField> (f: &DensePolynomial<F>, gamma: &F) -> DensePolynomial<F>{
    let mut coeffs = f.coeffs().to_vec();
    let mut curr_gamma = F::one();
    for i in 0..coeffs.len(){
        coeffs[i] *= curr_gamma;
        curr_gamma *= gamma;
    }    DensePolynomial::<F>::from_coefficients_vec(coeffs)
}
#[allow(non_snake_case)]
// return `x^n * f(x)` 
pub fn x_power_n_times_f<F: PrimeField>(f: &DensePolynomial<F>, N: usize) -> DensePolynomial<F>{
    let mut coeffs = vec![F::zero(); N];
    coeffs.append(&mut f.coeffs().to_vec());

    DensePolynomial::<F>::from_coefficients_vec(coeffs)
}
#[allow(non_snake_case)]
/// return hadamard product of f_1 and f_2
pub fn hadamard_product<F: PrimeField>(f_1: &DensePolynomial<F>, f_2: &DensePolynomial<F>) -> DensePolynomial<F>{
    let coeffs_1 = f_1.coeffs();
    let coeffs_2 = f_2.coeffs();

    let size = min(coeffs_1.len(), coeffs_2.len());
    let mut coeffs = vec![F::zero(); size];
    for i in 0..size{
        coeffs[i] = coeffs_1[i] * coeffs_2[i];
    }

    DensePolynomial::<F>::from_coefficients_vec(coeffs)
}

#[allow(non_snake_case)]
/// return 1 + x + x^2 + ... x^{size -1}
pub fn all_coeff_one_poly<F: PrimeField>(size: usize) -> DensePolynomial<F>{
    DensePolynomial::<F>::from_coefficients_vec(vec![F::one(); size])
}

/// return 1 + 2*x + 3*x^2 ... (size-1) x^(size-1)
#[allow(non_snake_case)]
pub fn id_poly<F: PrimeField>(size: usize) -> DensePolynomial<F> {
    let coeffs: Vec<F> = (0..size).map(|e| F::from(e as u64)).collect();
    DensePolynomial::from_coefficients_vec(coeffs)
}

/// applies the permutation sigma to polynomial f
#[allow(non_snake_case)]
pub fn compute_f_sigma<F: PrimeField>(f: &DensePolynomial<F>, sigma: &Vec<usize>) -> DensePolynomial<F>{
    let mut coeffs = vec![F::zero(); sigma.len()];
    let mut f_coeffs = f.coeffs().to_vec();
    while f_coeffs.len() < sigma.len(){
        f_coeffs.push(F::zero());
    }

    for i in 0..sigma.len(){
        coeffs[i] = f_coeffs[sigma[i]];
    }

    DensePolynomial::<F>::from_coefficients_vec(coeffs)
}

#[allow(non_snake_case)]
/// returns `x*f(x) mod x^N - 1`. Assumes that deg(f) < N
pub fn cyclic_right_shift<F: PrimeField>(f: &DensePolynomial<F>, N: usize) -> DensePolynomial<F>{
    let mut coeffs = f.coeffs().to_vec();
    let coeff_N = if coeffs.len() == N { coeffs[N-1] } else { F::zero() };
    let mut fin_coeffs = vec![coeff_N];
    fin_coeffs.append(&mut coeffs);
    if fin_coeffs.len() > N { fin_coeffs.pop(); };
    DensePolynomial::from_coefficients_vec(fin_coeffs)
}

/// add two points in affine form
pub fn add<E: Pairing>(p: &E::G1Affine, q: &E::G1Affine ) -> E::G1Affine {
    (p.into_group() + q.into_group()).into_affine()
}


/// subtract two points in affine form
pub fn sub<E: Pairing>(p: &E::G1Affine, q: &E::G1Affine ) -> E::G1Affine {
    (p.into_group() - q.into_group()).into_affine()
}

pub fn from<E: Pairing>(p: &Vec<E::G1Affine>) -> Vec<E::G1Prepared>
{
    p.iter().map(|q| E::G1Prepared::from(q)).collect()
}

/// converts field element to polynomial
pub fn as_poly<F: PrimeField>(x: F) -> DensePolynomial<F> {
    DensePolynomial::from_coefficients_vec(vec![x])
}

/// divides the polynomial by x, only when you can
pub fn divide_by_x<F: PrimeField>(f: &DensePolynomial<F>) -> DensePolynomial<F> {
    assert!(f.coeffs()[0] == F::zero());
    DensePolynomial::from_coefficients_slice(&f.coeffs()[1..])
}

/// parallel poly mul
pub fn parallel_poly_mul<F: PrimeField>(f: &[DensePolynomial<F>], g: &[DensePolynomial<F>]) -> Vec<DensePolynomial<F>> {
    f.par_iter().zip(g.par_iter()).map(|(f_i, g_i)| mul_ark_poly(f_i, g_i)).collect()
}