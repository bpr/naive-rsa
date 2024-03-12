use num_bigint::{BigInt, BigUint, ToBigInt, RandBigInt};
use num_traits::Zero;
use std::{cmp::PartialEq, ops::{Add, Mul, Neg}};

// A macro to create a BigInt from a string literal
#[macro_export]
macro_rules! bi {
    ($x: expr) => {
        $crate::bi!($x, 10)
    };
    ($x: expr, $base: literal) => {
        num_bigint::BigInt::parse_bytes($x.as_bytes(), $base).unwrap()
    }
}


pub struct PublicKey {
    n: BigInt,
    e: BigInt,
}

pub struct PrivateKey {
    d: BigInt,
}

pub fn is_even(n: BigInt) -> bool {
    n % 2 == BigInt::from(0)
}

pub fn factor_out_twos(n: BigInt) -> (usize, BigInt) {
    let mut s = 0;
    let mut d = n.clone();
    while is_even(d.clone()) {
        d /= 2;
        s += 1;
    }
    (s, d)
}

// is_probable_prime: determine if a number is probably prime using Miller-Rabin test
pub fn is_probable_prime(n: BigInt, num_rounds: usize) -> bool {
    // If n is even, it's not prime
    if is_even(n.clone()) {
        return false;
    }
    let (s, d) = factor_out_twos(n.clone() - 1);
    let n1: BigInt = n.clone() - 1;
    for _ in 0..num_rounds {
        let a = rand::thread_rng().gen_bigint_range(&BigInt::from(2), &n1);
        let mut x = a.modpow(&d, &n.clone());
        for _ in 0..s {
            let y = x.modpow(&BigInt::from(2), &n.clone());
            if y == BigInt::from(1) && x != BigInt::from(1) && x != n1 {
                return false;
            }
            x = y;
        }
        if x != BigInt::from(1) {
            return false;
        }
    }
    // If we haven't found a witness, then n is probably prime
    true
}

pub fn random_prime(ndigits: u32) -> BigInt {
    let mut rng = rand::thread_rng();
    let low = BigInt::from(10).pow(ndigits - 1 as u32);
    let high = low.clone().pow(2);
    let mut p = rng.gen_bigint_range(&low, &high);
    while !is_probable_prime(p.clone(), 100) {
        p = rng.gen_bigint_range(&low, &high);
    }
    p
}

pub fn gen_keys() -> (PublicKey, PrivateKey) {
    // Pick two large primes p and q
    let p: BigInt = random_prime(100u32);
    let q: BigInt = random_prime(100u32);
    // Compute n = pq
    let n: BigInt = &p * &q;
    // Compute (p-1)(q-1)
    let phi: BigInt = (&p - 1) * (&q - 1);
    // Find 'e' relatively prime to (p-1)(q-1)
    let e: BigInt = BigInt::from(65537);
    let d = mod_inverse(e, phi.clone());
    (PublicKey { 
        n: n.clone(), 
        e: BigInt::from(65537),
     },
     PrivateKey {
        d,
     })
}

pub fn encrypt(pub_key: &PublicKey, m: BigInt) -> BigInt {
    m.modpow(&pub_key.e, &pub_key.n)
}

pub fn decrypt(pub_key: &PublicKey, priv_key: &PrivateKey, c: BigInt) -> BigInt {
    c.modpow(&priv_key.d, &pub_key.n)
}

// extended gcd: https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm a helpful utility function
pub fn extended_gcd(a: BigInt, b: BigInt) -> (BigInt, (BigInt, BigInt), (BigInt, BigInt)) {
    let (mut old_r, mut r) = (a, b);
    let (mut old_s, mut s) = (BigInt::from(1), BigInt::from(0));
    let (mut old_t, mut t) = (BigInt::from(0), BigInt::from(1));
    
    while r != BigInt::zero() {
        let quotient = old_r.clone() / r.clone();
        let tmp = r.clone();
        (old_r, r) = (r, old_r - quotient.clone() * tmp);
        let tmp = s.clone();
        (old_s, s) = (s, old_s - quotient.clone() * tmp);
        let tmp = t.clone();
        (old_t, t) = (t, old_t - quotient.clone() * tmp);
    }
    (old_r, (old_s, old_t), (s, t))
}

pub fn mod_inverse(a: BigInt, m: BigInt) -> BigInt {
    let (gcd, (s, _), _) = extended_gcd(a.clone(), m.clone());
    if gcd != BigInt::from(1) {
        panic!("{} and {} are not coprime", a, m);
    }
    if s < BigInt::zero() {
        m + s
    } else {
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const P: &str = "fffffffffffffffffffffffffffffffffffffffffffffffffffffffeffffffffffffffffffffffffffffffffffffffffffffffffffffffff";

    #[test]
    fn primality_test_works_on_primes() {
        assert!(is_probable_prime(bi!("17"), 20));
        assert!(is_probable_prime(bi!(P, 16), 20));
    }

    #[test]
    fn primality_test_works_on_composites() {
        assert!(!is_probable_prime(bi!("355") * bi!("113"), 20));
    }

    #[test]
    fn encryption_and_decryption_work_on_u8() {
        let (pub_key, priv_key) = gen_keys();
        for i in 0..255 {
            let m = BigInt::from(i);
            let c = encrypt(&pub_key, m.clone());
            let m_prime = decrypt(&pub_key, &priv_key, c);
            assert_eq!(m, m_prime);
        }
    }
}
