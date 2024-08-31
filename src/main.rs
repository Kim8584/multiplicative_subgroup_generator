mod sanity_checks {
    use rand::Rng;
    // check that a number is a factor to another number
    pub fn is_factor(factor: i32, number: i32) -> bool {
        number % factor == 0
    }
    pub fn mod_exp(mut a: u64, mut s: u64, n: u64) -> u64 {
        let mut result = 1;
        a %= n;
        while s > 0 {
            if s % 2 == 1 {
                result = (result * a) % n;
            }
            s /= 2;
            a = (a * a) % n;
        }
        result
    }
    // check if a number is prime using miller rabin algo
    pub fn is_prime(n: u64, k: u64) -> bool {
        if n <= 1 || n == 4 {
            return false;
        }
        if n <= 3 {
            return true;
        }

        let mut r = 0;
        let mut s = n - 1;
        while s % 2 == 0 {
            r += 1;
            s /= 2;
        }

        let mut rng = rand::thread_rng();
        for _ in 0..k {
            let a: u64 = rng.gen_range(2..n - 2);
            let mut x = mod_exp(a, s, n);
            if x == 1 || x == n - 1 {
                continue;
            }
            let mut is_composite = true;
            for _ in 0..r - 1 {
                x = mod_exp(x, 2, n);
                if x == n - 1 {
                    is_composite = false;
                    break;
                }
            }
            if is_composite {
                return false;
            }
        }
        true
    }
}
mod primitive_root {
    use crate::sanity_checks::mod_exp;
    // find factors of k
    pub fn factors(k: u64) -> Vec<u64> {
        let mut factors = Vec::new();
        for i in 1..=k {
            if k % i == 0 {
                factors.push(i);
            }
        }
        factors
    }
}
mod multiplicative_subgruop {
    use crate::error::*;
    use crate::field::{generate_candidate, is_generator};
    use crate::primitive_root::factors;
    use crate::sanity_checks::mod_exp;
    use rand::Rng;
    use std::collections::HashSet;
    use std::iter::FromIterator;

    // generate the multiplicative subgroup of size n from field modulo p
    // this function returns the multplicative subgroup of size n from field modulo p
    // it first checks if p is prime
    // then it checks than n is a factor of p-1
    // p is not prime it returns an error and if n is not a factor of p-1 it returns an error
    // then it generates a candidate for the primitive root
    // then it checks if the candidate is a primitive root
    // if the candidate is a primitive root then it returns the multiplicative subgroup
    pub fn multiplicative_subgroup(p: u64, n: u64) -> Result<Vec<u64>, Box<dyn std::error::Error>> {
        if !crate::sanity_checks::is_prime(p, 5) {
            return Err(Box::new(NotPrimeError));
        }

        if (p - 1) % n != 0 {
            return Err(Box::new(NotFactorError));
        }
        // let mut rng = rand::thread_rng();
        let mut g = generate_candidate(p);
        while !is_generator(p, g) {
            g = generate_candidate(p);
        }
        let mut subgroup = HashSet::new();

        // we generate element in the subgroup by raising the generator to the power i((p-1)/n) mod p where i is in the range of 1 to n
        for i in 1..=n {
            subgroup.insert(mod_exp(g, i * ((p - 1) / n), p));
        }

        let mut subgroup = Vec::from_iter(subgroup);
        // rotate the list  until 1 is the first element in the list
        let index = subgroup.iter().position(|&x| x == 1).unwrap();
        subgroup.rotate_left(index);
        Ok(subgroup)
    }
}
// this mod is where i put error
mod error {
    // custom error if n is not a factor of p-1
    #[derive(Debug)]
    pub struct NotFactorError;
    impl std::error::Error for NotFactorError {
        fn description(&self) -> &str {
            "n is not a factor of p-1"
        }
    }
    impl std::fmt::Display for NotFactorError {
        fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            write!(f, "n is not a factor of p-1")
        }
    }
    //custom error if p is not prime
    #[derive(Debug)]
    pub struct NotPrimeError;
    impl std::error::Error for NotPrimeError {
        fn description(&self) -> &str {
            "p is not prime"
        }
    }
    impl std::fmt::Display for NotPrimeError {
        fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            write!(f, "p is not prime")
        }
    }
}
mod field {
    use crate::sanity_checks::mod_exp;
    use rand::Rng;
    // generate a random element from a field of modulo p the random element will be tested to be a valid primimitive root
    // to get a valid candidate the elements are in the range of 2 <= x <= p-2
    pub fn generate_candidate(p: u64) -> u64 {
        let mut rng = rand::thread_rng();
        rng.gen_range(2..p - 1)
    }
    // check if a number is a primitive root modulo p
    pub fn is_generator(p: u64, g: u64) -> bool {
        let mut factors = crate::primitive_root::factors(p - 1);
        // pop the last element of the factors since it is p-1
        factors.pop();
        // remove the first element in factors since it is 1
        factors.remove(0);
        // check if g is a primitive root modulo p using the factors
        // if g^((p-1)/f) mod p == 1 for all factors f of p-1 then g is not a primitive root modulo p
        for f in factors {
            if mod_exp(g, (p - 1) / f, p) == 1 {
                return false;
            }
        }
        true
    }
}

fn main() {
    println!("Hello, world!");
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::multiplicative_subgruop::multiplicative_subgroup;
    use crate::primitive_root::factors;
    use field::is_generator;

    #[test]
    fn test_factors() {
        assert_eq!(factors(12), vec![1, 2, 3, 4, 6, 12]);
        assert_eq!(factors(15), vec![1, 3, 5, 15]);
        assert_eq!(factors(17), vec![1, 17]);
        assert_eq!(factors(24), vec![1, 2, 3, 4, 6, 8, 12, 24]);
        assert_eq!(factors(6), vec![1, 2, 3, 6]);
    }
    // test is_generator function
    #[test]
    fn test_is_generator() {
        assert_eq!(is_generator(7, 3), true);
        assert_eq!(is_generator(11, 2), true);
        assert_eq!(is_generator(13, 2), true);
        assert_eq!(is_generator(17, 3), true);
        assert_eq!(is_generator(19, 2), true);
        assert_eq!(is_generator(23, 5), true);
        assert_eq!(is_generator(29, 2), true);
        assert_eq!(is_generator(31, 3), true);
        assert_eq!(is_generator(37, 2), true);
        assert_eq!(is_generator(41, 6), true);
        assert_eq!(is_generator(43, 3), true);
        assert_eq!(is_generator(47, 5), true);
        assert_eq!(is_generator(53, 2), true);
        assert_eq!(is_generator(59, 2), true);
        assert_eq!(is_generator(61, 2), true);
        assert_eq!(is_generator(67, 2), true);
        assert_eq!(is_generator(71, 7), true);
        assert_eq!(is_generator(73, 5), true);
        assert_eq!(is_generator(79, 3), true);
        assert_eq!(is_generator(83, 2), true);
        assert_eq!(is_generator(89, 3), true);
        assert_eq!(is_generator(97, 5), true);
        assert_eq!(is_generator(101, 2), true);
        assert_eq!(is_generator(103, 5), true);
        assert_eq!(is_generator(107, 2), true);
        assert_eq!(is_generator(109, 6), true);
        assert_eq!(is_generator(113, 3), true);
        assert_eq!(is_generator(127, 3), true);
        // assert_eq!(is_generator(337, 85), true);
    }
    // test miller rabin working correctly so test is prime
    #[test]
    fn test_is_prime() {
        assert_eq!(crate::sanity_checks::is_prime(7, 5), true);
        assert_eq!(crate::sanity_checks::is_prime(11, 5), true);
        assert_eq!(crate::sanity_checks::is_prime(13, 5), true);
        assert_eq!(crate::sanity_checks::is_prime(17, 5), true);
        assert_eq!(crate::sanity_checks::is_prime(19, 5), true);
        assert_eq!(crate::sanity_checks::is_prime(23, 5), true);
        assert_eq!(crate::sanity_checks::is_prime(29, 5), true);
        assert_eq!(crate::sanity_checks::is_prime(31, 5), true);
        assert_eq!(crate::sanity_checks::is_prime(37, 5), true);
        assert_eq!(crate::sanity_checks::is_prime(41, 5), true);
        assert_eq!(crate::sanity_checks::is_prime(43, 5), true);
        assert_eq!(crate::sanity_checks::is_prime(47, 5), true);
        assert_eq!(crate::sanity_checks::is_prime(53, 5), true);
        assert_eq!(crate::sanity_checks::is_prime(59, 5), true);
        assert_eq!(crate::sanity_checks::is_prime(61, 5), true);
        assert_eq!(crate::sanity_checks::is_prime(67, 5), true);
        assert_eq!(crate::sanity_checks::is_prime(71, 5), true);
        assert_eq!(crate::sanity_checks::is_prime(73, 5), true);
        assert_eq!(crate::sanity_checks::is_prime(79, 5), true);
        assert_eq!(crate::sanity_checks::is_prime(83, 5), true);
        assert_eq!(crate::sanity_checks::is_prime(89, 5), true);
        assert_eq!(crate::sanity_checks::is_prime(97, 5), true);
        assert_eq!(crate::sanity_checks::is_prime(101, 5), true);
        assert_eq!(crate::sanity_checks::is_prime(103, 5), true);
    }
    // test the multiplicative subgroup
    // asserts that for functions with n not a foctor of p - 1 returns error
    #[test]
    fn test_multiplicative_subgroup() {
        assert_eq!(multiplicative_subgroup(7, 3).unwrap(), vec![1, 2, 4]);
        // assert_eq!(multiplicative_subgroup(11, 5).unwrap(), vec![1, 3, 4, 5, 9]);
    }
}
