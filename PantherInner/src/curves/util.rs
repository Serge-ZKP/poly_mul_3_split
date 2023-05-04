
/// Since BigUint allows only conversion from u32 limbs, this helper function
/// is implemented to convert from u64 limbs to u32 limbs.
pub fn u64_to_u32(slice: &[u64],arr: &mut [u32])->(){
    for i in 0..slice.len(){
        arr[2*i]=slice[i] as u32;
        arr[2*i+1]=(slice[i]>>32) as u32;
    }
}

#[test]
fn test_u64_to_u32(){
    use num_bigint::BigUint;
    use ark_std::UniformRand;
    let mut rng=ark_std::test_rng();
    //let num= "988340759509845985983392897970799781".parse::<BigUint>().unwrap(); //example to parse string to BigUint
    let num_as_big_int:ark_ff::BigInt<10>= ark_ff::BigInt::rand(&mut rng);
    let num_as_big_uint=BigUint::from(num_as_big_int);
    let u64rep = num_as_big_uint.to_u64_digits();
    let mut u32rep=vec!(0 as u32;2*u64rep.len());
    u64_to_u32(&u64rep[..],&mut u32rep[..]);
    assert_eq!(BigUint::from_slice(&u32rep[..]),num_as_big_uint)
}