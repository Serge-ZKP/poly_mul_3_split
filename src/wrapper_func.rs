
#![allow(non_snake_case)]
#![allow(dead_code)]

use std::ffi::CString;
use std::sync::Mutex;
use ark_ff::BigInteger;
use ark_ff::PrimeField;
use ark_poly::DenseUVPolynomial;
use ark_poly::Polynomial;
use ark_poly::univariate::DenseOrSparsePolynomial;
use ark_poly::{self, univariate::DensePolynomial};

use rug::Integer;
use rug::integer::Order;
use wrapper_bindings::*;
use rayon::prelude::*;


/// Takes a field element and initiates the NTL modulus.
pub fn ark_NTL_modulus_init<F: ark_ff::PrimeField>() {
    let modulus = Integer::from_digits(F::characteristic(), Order::Lsf); 
    let modulus_string =  CString::new(modulus.to_string()).unwrap();
    unsafe{
        let modulus = ZZ_set_string(modulus_string.as_ptr() as *const i8);
        init_modulus(modulus);
        ZZ_destroy(modulus);
        
    }
}

/// Converts arkworks Field element to an NTL ZZ_p.
pub fn ark_to_NTL_ZZ_p<F: ark_ff::PrimeField>(x: F, ntl_x: &Send_ZZ_p) {
    let x_as_integer = Integer::from_digits(x.into_bigint().to_bytes_le().as_slice(), Order::Lsf);
    let x_cstr = CString::new(x_as_integer.to_string()).unwrap();

    unsafe {
        ZZ_p_set_string(x_cstr.as_ptr(), ntl_x.0);
    }
}

/// set coefficient of polynomial at given index 
pub fn set_coeff(poly: &Send_ZZ_pX, index: usize, coeff: &Send_ZZ_p){
    let index_as_i64 = index as i64;
    unsafe{
        ZZ_pX_SetCoeff_fp(poly.0, index_as_i64, coeff.0);
        //ZZ_p_destroy(coeff.0);
    }
}

/// returns coefficient of the polynomial at given index
pub fn get_coeff(poly: &Send_ZZ_pX, index: i64) -> *const ZZ_p{
    unsafe{
        ZZ_pX_coeff(poly.0, index)
    }
}
/// Converts an arkworks polynomial to an NTL ZZ_pX.
pub fn ark_to_NTL_ZZ_pX<F: ark_ff::PrimeField>(ark_poly: &DensePolynomial<F>, ntl_poly: &Send_ZZ_pX) {

    let coeff_list: Vec<F> = ark_poly.coeffs().to_vec();

    unsafe {
        let mutex_for_ntl_poly = Mutex::new(ntl_poly);
        coeff_list.par_iter().enumerate().for_each(|(i, e) | {
            ark_NTL_modulus_init::<F>();
            let coeff  = Send_ZZ_p(ZZ_p_zero());
            ark_to_NTL_ZZ_p(*e, &coeff);
            {
                let locked_mutex_for_ntl_poly = mutex_for_ntl_poly.lock().unwrap();
                set_coeff(*locked_mutex_for_ntl_poly, i, &coeff);
                ZZ_p_destroy(coeff.0);
            };
        });
    }
}


/// Converts an NTL ZZ_p to an arkworks field element.
pub fn NTL_ZZ_p_to_ark<F: ark_ff::PrimeField>(ntl_x: *const ZZ_p) -> F {
    let len = Box::into_raw(Box::new(0 as i64));
    unsafe{
        let c_ptr = ZZ_p_to_string(ntl_x, len);
        let x_slice = std::slice::from_raw_parts(c_ptr, *len as usize);
        let ark_elem = F::from_le_bytes_mod_order(x_slice);
        destroy_unsigned_char(c_ptr);
        ark_elem
    }
}

/// Converts an NTL ZZ_pX to an arkworks polynomial.
pub fn NTL_ZZ_pX_to_ark<F: ark_ff::PrimeField>(ntl_poly: Send_ZZ_pX) -> DensePolynomial<F> 
{
    unsafe {
        let deg = ZZ_pX_deg(ntl_poly.0);
        let mutex_for_ntl_poly = Mutex::new(&ntl_poly);
        let coeffs = (0..deg+1).collect::<Vec<i64>>().par_iter().map(|i| {
            let coeff = {
                let locked_mutex_for_ntl_poly = mutex_for_ntl_poly.lock().unwrap();
                get_coeff(&*locked_mutex_for_ntl_poly, *i)
            };
            NTL_ZZ_p_to_ark(coeff)
        }).collect();
        ZZ_pX_destroy(ntl_poly.0);
        DensePolynomial::from_coefficients_vec(coeffs)
    }
    
}

/// multiply two ark polynomials. Uses highly optimized NTL functions under the hood
pub fn mul_ark_poly<F: PrimeField>(f: &DensePolynomial<F>, g: &DensePolynomial<F>) -> DensePolynomial<F>{
    unsafe {
        ark_NTL_modulus_init::<F>();
        let x = Send_ZZ_pX(ZZ_pX_zero());
        let y = Send_ZZ_pX(ZZ_pX_zero());
        let z = Send_ZZ_pX(ZZ_pX_zero());

        ark_to_NTL_ZZ_pX(f, &x);
        ark_to_NTL_ZZ_pX(g, &y);

        ZZ_pX_mul(z.0, x.0, y.0);

        ZZ_pX_destroy(x.0);
        ZZ_pX_destroy(y.0);

        let poly = NTL_ZZ_pX_to_ark(z);
        poly
    }
}

/// sum of products of several polynomials
pub fn sum_of_prods<F: PrimeField>(vec_fg: Vec<(&DensePolynomial<F>, &DensePolynomial<F>)>) -> DensePolynomial<F>{
    unsafe{
        ark_NTL_modulus_init::<F>();
        let sum = vec_fg.par_iter().map(|(f, g)| {
            
            let x = Send_ZZ_pX(ZZ_pX_zero());
            let y = Send_ZZ_pX(ZZ_pX_zero());
            let z = Send_ZZ_pX(ZZ_pX_zero());
    
            ark_to_NTL_ZZ_pX(f, &x);
            ark_to_NTL_ZZ_pX(g, &y);
    
            ZZ_pX_mul(z.0, x.0, y.0);
            z
        }).reduce(|| { Send_ZZ_pX(ZZ_pX_zero()) }, |acc, h| {  ZZ_pX_add(acc.0, h.0, acc.0); acc});
        NTL_ZZ_pX_to_ark(sum)
    }
}
/// Divide f by linear polynomial g and return quotient and remainder
pub fn div_by_linear_poly<F: PrimeField>(f: &DensePolynomial<F>, g: &DensePolynomial<F>) -> (DensePolynomial<F>, F){
    assert!(g.degree() <= 1);
    let (q, r) = DenseOrSparsePolynomial::from(f).divide_with_q_and_r(&DenseOrSparsePolynomial::from(g)).unwrap();
    (q, if r.coeffs().len() > 0 {r.coeffs()[0]} else {F::zero()} )
}

/// used for debugging purposes.
pub fn print_poly<F: PrimeField>(f: &DensePolynomial<F>){
    unsafe{
        let x = Send_ZZ_p(ZZ_p_zero());
        for e in f.coeffs(){
            ark_to_NTL_ZZ_p(*e, &x);
            ZZ_p_print(x.0);
        }
    }
   
}
#[cfg(test)]
mod wrapper_func_tests{
    use crate::{wrapper_func::*, util::rand_poly};
    use ark_ff::Zero;
    use ark_std::UniformRand;
    //use ark_bls12_381::Bls12_381;
    use panther_inner::curves::Panther12_598;
    use ark_ec::pairing::Pairing;
    type F = <Panther12_598 as Pairing>::ScalarField;
    use ark_std::rand::rngs::StdRng;
    use rand_core::SeedableRng;
    use std::{time::Instant};
    macro_rules! init_test  {
        ($rng: ident, $deg: ident, $num: expr) => {
            let $rng = &mut StdRng::from_entropy();
            let $deg = $num as usize;
            ark_NTL_modulus_init::<F>();
        };
    }
    
    #[test]
    fn test_conv_between_ark_to_ntl_ZZ_p(){
        init_test!(rng, deg, 10000);
        unsafe{
            let now = Instant::now();
            for _ in 0..deg{
                let p = F::rand(rng);
                let x = Send_ZZ_p(ZZ_p_zero());
                ark_to_NTL_ZZ_p(p, &x);
                let q: F= NTL_ZZ_p_to_ark(x.0);
                assert_eq!(p, q);
            }
            let elapsed = now.elapsed();
            println!("Elapsed: {:.2?}", elapsed);
        }
    }
    
    #[test]
    fn time_ark_to_ntl_poly(){
        init_test!(rng, deg, 1000000);
        unsafe{
            let mut coeffs = vec![];
            for _ in 0..deg{
                let p = F::rand(rng);
                coeffs.push(p);
            }
            let poly = DensePolynomial::<F>::from_coefficients_vec(coeffs);
            let xpoly = Send_ZZ_pX(ZZ_pX_zero());
            let now = Instant::now();
            ark_to_NTL_ZZ_pX(&poly, &xpoly);
            let elapsed = now.elapsed();
            println!("Elapsed: {:.2?}", elapsed);
        }
    }
    
    #[test]
    fn time_ntl_to_ark_poly(){
        init_test!(rng, deg, 1000000);
        unsafe{
            let mut coeffs = vec![];
            for _ in 0..deg{
                let p = F::rand(rng);
                coeffs.push(p);
            }
            let f = DensePolynomial::<F>::from_coefficients_vec(coeffs);
            let xpoly = Send_ZZ_pX(ZZ_pX_zero());
            ark_to_NTL_ZZ_pX(&f, &xpoly);
            let now = Instant::now();
            let _g = NTL_ZZ_pX_to_ark::<F>(xpoly);
            let elapsed = now.elapsed();
            println!("Elapsed: {:.2?}", elapsed);
        }
    }
    
    #[test]
    fn test_conv_between_ark_and_ntl_ZZ_pX(){
        init_test!(rng, deg, 1000000);
        unsafe{
            let mut coeffs = vec![];
            for _ in 0..deg{
                let p = F::rand(rng);
                coeffs.push(p);
            }
            let f = DensePolynomial::<F>::from_coefficients_vec(coeffs);
            let now = Instant::now();
            let xpoly = Send_ZZ_pX(ZZ_pX_zero());
            ark_to_NTL_ZZ_pX(&f, &xpoly);
            let g = NTL_ZZ_pX_to_ark::<F>(xpoly);
            assert_eq!(f, g);
            let elapsed = now.elapsed();
            println!("Elapsed: {:.2?}", elapsed);
        }
    }
    #[test]
    fn test_mul(){
        init_test!(rng, deg, 100);
        let (mut coeffs1, mut coeffs2) = (vec![], vec![]);
        for _ in 0..deg{
            coeffs1.push(F::rand(rng));
            coeffs2.push(F::rand(rng));
        }
    
        let f = DensePolynomial::<F>::from_coefficients_vec(coeffs1);
        let g = DensePolynomial::<F>::from_coefficients_vec(coeffs2);
        let actual = f.naive_mul(&g);
        let found = mul_ark_poly(&f, &g);
        assert_eq!(actual, found);
    
    }


    #[test]
    fn time_mul(){
        init_test!(rng, deg, 3000000);
        let (mut coeffs1, mut coeffs2) = (vec![], vec![]);
        for _ in 0..deg{
            coeffs1.push(F::rand(rng));
            coeffs2.push(F::rand(rng));
        }
        
        let f = DensePolynomial::<F>::from_coefficients_vec(coeffs1);
        let g = DensePolynomial::<F>::from_coefficients_vec(coeffs2);
    
        let now = Instant::now();
        let _h = mul_ark_poly(&f, &g);
        let elapsed = now.elapsed();
        println!("Elapsed: {:.2?}", elapsed);
    }

    /*
    Split both into 3 parts.

    f(X) = f_1(X) + f_2(X) X^n + f_3(X) X^{2n}

    h(X) = h_1(X) + h_2(X) X^n + h_3(X) X^{2n}

    Now compute the nine polyonmials f_i(X) * h_j(X)

    Then obtain f(X) * h(X) as the sum \sum_{i,j} [f_i(X) * h_j(X)] X^{i+j-2}
    sum_{i,j} [f_i(X) * h_j(X)] X^{n(i+j-2)}

    of 9 polynomials
    */
    #[test]
    fn time_mul_new() {
        //init_test!(rng, deg, 3000000);
        init_test!(rng, deg, 1000000);

        let mut f_coeff = vec![vec![]];
        let mut g_coeff = vec![vec![]];

        for _ in 0..3 {
            let (mut coeffs1, mut coeffs2) = (vec![], vec![]);
            for _ in 0..deg {
                coeffs1.push(F::rand(rng));
                coeffs2.push(F::rand(rng));
            }
            f_coeff.push(coeffs1);
            g_coeff.push(coeffs2);
        }

        let mut f = vec![];
        let mut g = vec![];

        for i in 0..3 {
                f.push( DensePolynomial::<F>::from_coefficients_vec(f_coeff[i].clone()));
                g.push( DensePolynomial::<F>::from_coefficients_vec(g_coeff[i].clone()));
        }

        let vec_fg = f.iter().zip(g.iter()).collect();

        let now = Instant::now();
        let h = sum_of_prods::<F>(vec_fg);
        let elapsed = now.elapsed();
        println!("Elapsed: {:.2?}", elapsed);


    }


    #[test]
    fn test_div(){
        init_test!(rng, deg, 1000000);
        let mut coeffs = vec![];
        for _ in 0..deg{
            let p = F::rand(rng);
            coeffs.push(p);
        }

        
        let f = DensePolynomial::from_coefficients_vec(coeffs);
        let lin = DensePolynomial::from_coefficients_vec(vec![F::rand(rng), F::rand(rng)]);

        let now = Instant::now();
        let (q, r) = div_by_linear_poly(&f, &lin);
        let elapsed = now.elapsed();
        println!("Elapsed: {:.2?}", elapsed);

        assert_eq!(q.naive_mul(&lin) + DensePolynomial::from_coefficients_vec(vec![r]), f);
    }

    #[test]
    fn test_sum_of_prods(){
        let n = 6;
        init_test!(_rng, deg, 100);
        let mut f = vec![];
        let mut g = vec![];
        for _ in 0..n{
            f.push(rand_poly(deg));
            g.push(rand_poly(deg));
        }
        let vec_fg = f.iter().zip(g.iter()).collect::<Vec<(&DensePolynomial<F>, &DensePolynomial<F>)>>();

        let ans = sum_of_prods::<F>(vec_fg.clone());
        let actual_ans = vec_fg.iter().map(|(f, g)| f.naive_mul(g)).fold(DensePolynomial::zero(), |acc, h| acc + h);
        assert_eq!(ans, actual_ans);
    }

    #[test]
    fn time_sum_of_prods(){
        let n = 6;
        init_test!(_rng, deg, 1000000);
        let mut f = vec![];
        let mut g = vec![];
        for _ in 0..n{
            f.push(rand_poly(deg));
            g.push(rand_poly(deg));
        }
        let vec_fg = f.iter().zip(g.iter()).collect();

        let now = Instant::now();
        let _ans = sum_of_prods::<F>(vec_fg);
        let elapsed = now.elapsed();
        println!("Elapsed: {:.2?}", elapsed);
    }
    
}
