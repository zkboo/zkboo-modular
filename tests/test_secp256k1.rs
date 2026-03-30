// SPDX-License-Identifier: LGPL-3.0-or-later

mod common;

use crate::common::rand_words::test_vec;
use core::array;
use dashu_int::UBig;
use zkboo::{
    executor::{ExecutionBackend, OwnedFlexibleWordPool},
    word::CompositeWord,
};
use zkboo_modular::montgomery::{MontgomeryMod, MontgomeryWord, MontgomeryWordRef};

const NUM_SAMPLES: usize = 100;

/// Modulus for the secp256k1 curve, implemented using two u128 limbs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SECP256K1;

impl MontgomeryMod<u128, 2> for SECP256K1 {
    #[inline]
    fn n(&self) -> CompositeWord<u128, 2> {
        return CompositeWord::from_be_words([
            0xffffffffffffffffffffffffffffffff,
            0xfffffffffffffffffffffffefffffc2f,
        ]);
    }

    #[inline]
    fn inv_exp(&self) -> Option<CompositeWord<u128, 2>> {
        return Some(CompositeWord::from_be_words([
            0xffffffffffffffffffffffffffffffff,
            0xfffffffffffffffffffffffefffffc2d,
        ]));
    }

    #[inline]
    fn rr_mod_n(&self) -> CompositeWord<u128, 2> {
        return CompositeWord::from_be_words([0x0, 0x1000007a2000e90a1]);
    }

    #[inline]
    fn n_neg_inv(&self) -> CompositeWord<u128, 2> {
        return CompositeWord::from_be_words([
            0xc9bd1905155383999c46c2c295f2b761,
            0xbcb223fedc24a059d838091dd2253531,
        ]);
    }
}

macro_rules! test_unop {
    ($func: ident, |$in_: ident, $m: ident| $body: block) => {
        ::paste::paste! {
            #[test]
            fn [<test_ $func>](){
                {
                    let seed = 0u64;
                    let $m = SECP256K1;
                    let samples_in_ = test_vec::<u128, 2, CompositeWord<u128, 2>>(NUM_SAMPLES, seed+1);
                    for $in_ in samples_in_ {
                        $body
                    }
                }
            }
        }
    };
}

macro_rules! test_binop {
    ($func: ident, |$inl: ident, $inr: ident, $m: ident| $body: block) => {
        ::paste::paste! {
            #[test]
            fn [<test_ $func>](){
                {
                    let seed = 0u64;
                    let $m = SECP256K1;
                    let samples_inl = test_vec::<u128, 2, CompositeWord<u128, 2>>(NUM_SAMPLES, seed+1);
                    let samples_inr = test_vec::<u128, 2, CompositeWord<u128, 2>>(NUM_SAMPLES, seed+2);
                    for ($inl, $inr) in samples_inl.into_iter().zip(samples_inr.into_iter()) {
                        $body
                    }
                }
            }
        }
    };
}

type WP = OwnedFlexibleWordPool<usize>;

fn to_ubig(w: CompositeWord<u128, 2>) -> UBig {
    let [lo, hi] = w.to_le_words();
    return UBig::from(hi) << 128 | UBig::from(lo);
}

fn from_ubig(x: UBig) -> CompositeWord<u128, 2> {
    let lo = (&x & UBig::from(u128::MAX)).try_into().unwrap();
    let hi = (&x >> 128).try_into().unwrap();
    return CompositeWord::from_le_words([lo, hi]);
}

test_unop!(const_reduce, |in_, m| {
    let out = {
        let in_ = MontgomeryWord::new(in_, m);
        let out = in_;
        out.value()
    };
    let reference_out = {
        let in_ = to_ubig(in_);
        let n = to_ubig(m.n());
        let out = in_ % n;
        from_ubig(out)
    };
    assert_eq!(out, reference_out);
});

test_unop!(const_neg, |in_, m| {
    let out = {
        let in_ = MontgomeryWord::new(in_, m);
        let out = -in_;
        out.value()
    };
    let reference_out = {
        let in_ = to_ubig(in_);
        let n = to_ubig(m.n());
        let in_ = in_ % n.clone();
        let out = (n.clone() - in_) % n;
        from_ubig(out)
    };
    assert_eq!(out, reference_out);
});

test_binop!(const_add, |inl, inr, m| {
    let out = {
        let inl = MontgomeryWord::new(inl, m);
        let inr = MontgomeryWord::new(inr, m);
        let out = inl + inr;
        out.value()
    };
    let reference_out = {
        let inl = to_ubig(inl);
        let inr = to_ubig(inr);
        let n = to_ubig(m.n());
        let out = (inl + inr) % n;
        from_ubig(out)
    };
    assert_eq!(out, reference_out);
});

test_binop!(const_sub, |inl, inr, m| {
    let out = {
        let inl = MontgomeryWord::new(inl, m);
        let inr = MontgomeryWord::new(inr, m);
        let out = inl - inr;
        out.value()
    };
    let reference_out = {
        let inl = to_ubig(inl);
        let inr = to_ubig(inr);
        let n = to_ubig(m.n());
        let inl = inl % n.clone();
        let inr = inr % n.clone();
        let out = if inl >= inr { inl - inr } else { inl + n - inr };
        from_ubig(out)
    };
    assert_eq!(out, reference_out);
});

test_binop!(const_mul, |inl, inr, m| {
    let out = {
        let inl = MontgomeryWord::new(inl, m);
        let inr = MontgomeryWord::new(inr, m);
        let out = inl * inr;
        out.value()
    };
    let reference_out = {
        let inl = to_ubig(inl);
        let inr = to_ubig(inr);
        let n = to_ubig(m.n());
        let out = (inl * inr) % n;
        from_ubig(out)
    };
    assert_eq!(out, reference_out);
});

test_unop!(const_inv, |in_, m| {
    let out = {
        if in_.is_zero() {
            in_
        } else {
            let in_ = MontgomeryWord::new(in_, m);
            let out = in_ * in_.inv();
            out.value()
        }
    };
    let reference_out = if in_.is_zero() {
        in_
    } else {
        CompositeWord::<u128, 2>::ONE
    };
    assert_eq!(out, reference_out);
});

test_unop!(reduce, |in_, m| {
    let exec_out = {
        let executor = ExecutionBackend::<WP>::new().into_executor();
        let in_ = MontgomeryWordRef::new(executor.input(in_), m);
        let out = in_;
        executor.output(out.value());
        let outputs = executor.finalize().u128;
        CompositeWord::from_le_words(array::from_fn(|i| outputs[i]))
    };
    let reference_out = {
        let in_ = MontgomeryWord::new(in_, m);
        let out = in_;
        out.value().into()
    };
    assert_eq!(exec_out, reference_out);
});

test_unop!(neg, |in_, m| {
    let exec_out = {
        let executor = ExecutionBackend::<WP>::new().into_executor();
        let in_ = MontgomeryWordRef::new(executor.input(in_), m);
        let out = -in_;
        executor.output(out.value());
        let outputs = executor.finalize().u128;
        CompositeWord::from_le_words(array::from_fn(|i| outputs[i]))
    };
    let reference_out = {
        let in_ = MontgomeryWord::new(in_, m);
        let out = -in_;
        out.value().into()
    };
    assert_eq!(exec_out, reference_out);
});

test_binop!(add, |inl, inr, m| {
    let exec_out = {
        let executor = ExecutionBackend::<WP>::new().into_executor();
        let inl = MontgomeryWordRef::new(executor.input(inl), m);
        let inr = MontgomeryWordRef::new(executor.input(inr), m);
        let out = inl + inr;
        executor.output(out.value());
        let outputs = executor.finalize().u128;
        CompositeWord::from_le_words(array::from_fn(|i| outputs[i]))
    };
    let reference_out = {
        let inl = MontgomeryWord::new(inl, m);
        let inr = MontgomeryWord::new(inr, m);
        let out = inl + inr;
        out.value().into()
    };
    assert_eq!(exec_out, reference_out);
});

test_binop!(sub, |inl, inr, m| {
    let exec_out = {
        let executor = ExecutionBackend::<WP>::new().into_executor();
        let inl = MontgomeryWordRef::new(executor.input(inl), m);
        let inr = MontgomeryWordRef::new(executor.input(inr), m);
        let out = inl - inr;
        executor.output(out.value());
        let outputs = executor.finalize().u128;
        CompositeWord::from_le_words(array::from_fn(|i| outputs[i]))
    };
    let reference_out = {
        let inl = MontgomeryWord::new(inl, m);
        let inr = MontgomeryWord::new(inr, m);
        let out = inl - inr;
        out.value().into()
    };
    assert_eq!(exec_out, reference_out);
});

test_binop!(mul, |inl, inr, m| {
    let exec_out = {
        let executor = ExecutionBackend::<WP>::new().into_executor();
        let inl = MontgomeryWordRef::new(executor.input(inl), m);
        let inr = MontgomeryWordRef::new(executor.input(inr), m);
        let out = inl * inr;
        executor.output(out.value());
        let outputs = executor.finalize().u128;
        CompositeWord::from_le_words(array::from_fn(|i| outputs[i]))
    };
    let reference_out = {
        let inl = MontgomeryWord::new(inl, m);
        let inr = MontgomeryWord::new(inr, m);
        let out = inl * inr;
        out.value().into()
    };
    assert_eq!(exec_out, reference_out);
});

test_unop!(inv, |in_, m| {
    let exec_out = {
        let executor = ExecutionBackend::<WP>::new().into_executor();
        let in_ = MontgomeryWordRef::new(executor.input(in_), m);
        let out = in_.inv();
        executor.output(out.value());
        let outputs = executor.finalize().u128;
        CompositeWord::from_le_words(array::from_fn(|i| outputs[i]))
    };
    let reference_out = {
        let in_ = MontgomeryWord::new(in_, m);
        let out = in_.inv();
        out.value().into()
    };
    assert_eq!(exec_out, reference_out);
});
