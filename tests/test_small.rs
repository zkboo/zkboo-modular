// SPDX-License-Identifier: LGPL-3.0-or-later

mod common;

use crate::common::{rand_words::test_vec, utils::modinv_i128};
use zkboo::{
    executor::{ExecutionBackend, OwnedFlexibleWordPool},
    word::{CompositeWord, Word},
};
use zkboo_modular::montgomery::{MontgomeryMod, MontgomeryWord, MontgomeryWordRef};

const NUM_SAMPLES: usize = 100;

/// Modulus for the secp256k1 curve.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SmallMod<W: Word> {
    n: W,
}

impl<W: Word> SmallMod<W> {
    pub fn new(n: W) -> Self {
        assert!(
            W::WIDTH <= 64,
            "SmallMod only supports word sizes up to 64 bits."
        );
        assert!(n & W::ONE == W::ONE, "SmallMod only supports odd moduli.");
        return Self { n };
    }
}

impl<W: Word> MontgomeryMod<W, 1> for SmallMod<W> {
    #[inline]
    fn n(&self) -> CompositeWord<W, 1> {
        return self.n.into();
    }

    fn inv_exp(&self) -> Option<CompositeWord<W, 1>> {
        return None;
    }

    fn rr_mod_n(&self) -> CompositeWord<W, 1> {
        let n: u128 = self.n().into().cast();
        let r = 1u128 << W::WIDTH;
        let r_mod_n = r % n;
        return W::cast_from((r_mod_n * r_mod_n) % n).into();
    }

    #[inline]
    fn n_neg_inv(&self) -> CompositeWord<W, 1> {
        let n: u128 = self.n().into().cast();
        let r = 1u128 << W::WIDTH;
        let n_inv = modinv_i128(n as i128, r as i128).unwrap();
        assert!(n_inv > 0);
        assert!(n_inv < r as i128);
        let n_inv = n_inv as u128;
        assert!((n * n_inv) % r == 1, "n_inv is not correct");
        return W::cast_from(r.wrapping_sub(n_inv)).into();
    }
}

macro_rules! _on_all_words {
    ($m:ident ! ( $($args:tt)* )) => {
        $m!(
            [u8, u16, u32, u64],
            $($args)*
        );
    };
}

macro_rules! _test_unop {
    ([$($W_: ty),* $(,)?], $func: ident, $W: ident, |$in_: ident, $m: ident| $body: block) => {
        ::paste::paste!{
            #[test]
            fn [<test_ $func>](){
                $(
                    {
                        type $W = $W_;
                        assert!($W::MAX.trailing_ones() <= 64);
                        let seed = 0u64;
                        let mut samples_n = test_vec::<$W, 1, $W>(15, seed);
                        samples_n[0..5].copy_from_slice(&[3, 7, 21, ($W::MAX-1)/2, $W::MAX]);
                        for n in samples_n {
                            let $m = SmallMod::new(n | 1);
                            let samples_in_ = test_vec::<$W, 1, $W>(NUM_SAMPLES, seed+1);
                            for $in_ in samples_in_ {
                                $body
                            }
                        }
                    }
                )*
            }
        }
    };
}

macro_rules! test_unop {
    ($func: ident, $W: ident, |$in_: ident, $m: ident| $body: block) => {
        _on_all_words!(_test_unop!($func, $W, |$in_, $m| $body));
    };
}

macro_rules! _test_binop {
    ([$($W_: ty),* $(,)?], $func: ident, $W: ident, |$inl: ident, $inr: ident, $m: ident| $body: block) => {
        ::paste::paste!{
            #[test]
            fn [<test_ $func>](){
                $(
                    {
                        type $W = $W_;
                        assert!($W::MAX.trailing_ones() <= 64);
                        let seed = 0u64;
                        let mut samples_n = test_vec::<$W, 1, $W>(15, seed);
                        samples_n[0..5].copy_from_slice(&[3, 7, 21, ($W::MAX-1)/2, $W::MAX]);
                        for n in samples_n {
                            let $m = SmallMod::new(n | 1);
                            let samples_inl = test_vec::<$W, 1, $W>(NUM_SAMPLES, seed+1);
                            let samples_inr = test_vec::<$W, 1, $W>(NUM_SAMPLES, seed+2);
                            for ($inl, $inr) in samples_inl.into_iter().zip(samples_inr.into_iter()) {
                                $body
                            }
                        }
                    }
                )*
            }
        }
    };
}

macro_rules! test_binop {
    ($func: ident, $W: ident, |$inl: ident, $inr: ident, $m: ident| $body: block) => {
        _on_all_words!(_test_binop!($func, $W, |$inl, $inr, $m| $body));
    };
}

test_unop!(const_reduce, W, |in_, m| {
    assert_eq!(MontgomeryWord::new(in_, m).value().into(), {
        let n = m.n().into();
        in_ % n
    });
});

test_unop!(const_neg, W, |in_, m| {
    assert_eq!((-MontgomeryWord::new(in_, m)).value().into(), {
        let n = m.n().into();
        let in_ = in_ % n;
        if in_ == 0 { in_ } else { n - in_ }
    });
});

test_binop!(const_add, W, |inl, inr, m| {
    assert_eq!(
        {
            let inl = MontgomeryWord::new(inl, m);
            let inr = MontgomeryWord::new(inr, m);
            (inl + inr).value().into()
        },
        {
            let n = m.n().into();
            ((inl as u128 + inr as u128) % (n as u128)) as W
        }
    );
});

test_binop!(const_sub, W, |inl, inr, m| {
    assert_eq!(
        {
            let inl = MontgomeryWord::new(inl, m);
            let inr = MontgomeryWord::new(inr, m);
            (inl - inr).value().into()
        },
        {
            let n = m.n().into();
            let inl = inl % n;
            let inr = inr % n;
            if inl >= inr { inl - inr } else { n - inr + inl }
        }
    );
});

test_binop!(const_mul, W, |inl, inr, m| {
    assert_eq!(
        {
            let inl = MontgomeryWord::new(inl, m);
            let inr = MontgomeryWord::new(inr, m);
            (inl * inr).value().into()
        },
        {
            let n = m.n().into();
            ((inl as u128 * inr as u128) % (n as u128)) as W
        }
    );
});

type WP = OwnedFlexibleWordPool<usize>;

test_unop!(reduce, W, |in_, m| {
    let exec_out: W = {
        let executor = ExecutionBackend::<WP>::new().into_executor();
        let in_ = MontgomeryWordRef::new(executor.input(in_), m);
        let out = in_;
        executor.output(out.value());
        executor.finalize().as_vec()[0]
    };
    let reference_out = {
        let in_ = MontgomeryWord::new(in_, m);
        let out = in_;
        out.value().into()
    };
    assert_eq!(exec_out, reference_out);
});

test_unop!(neg, W, |in_, m| {
    let exec_out: W = {
        let executor = ExecutionBackend::<WP>::new().into_executor();
        let in_ = MontgomeryWordRef::new(executor.input(in_), m);
        let out = -in_;
        executor.output(out.value());
        executor.finalize().as_vec()[0]
    };
    let reference_out = {
        let in_ = MontgomeryWord::new(in_, m);
        let out = -in_;
        out.value().into()
    };
    assert_eq!(exec_out, reference_out);
});

test_binop!(add, W, |inl, inr, m| {
    let exec_out: W = {
        let executor = ExecutionBackend::<WP>::new().into_executor();
        let inl = MontgomeryWordRef::new(executor.input(inl), m);
        let inr = MontgomeryWordRef::new(executor.input(inr), m);
        let out = inl + inr;
        executor.output(out.value());
        executor.finalize().as_vec()[0]
    };
    let reference_out = {
        let inl = MontgomeryWord::new(inl, m);
        let inr = MontgomeryWord::new(inr, m);
        let out = inl + inr;
        out.value().into()
    };
    assert_eq!(exec_out, reference_out);
});

test_binop!(sub, W, |inl, inr, m| {
    let exec_out: W = {
        let executor = ExecutionBackend::<WP>::new().into_executor();
        let inl = MontgomeryWordRef::new(executor.input(inl), m);
        let inr = MontgomeryWordRef::new(executor.input(inr), m);
        let out = inl - inr;
        executor.output(out.value());
        executor.finalize().as_vec()[0]
    };
    let reference_out = {
        let inl = MontgomeryWord::new(inl, m);
        let inr = MontgomeryWord::new(inr, m);
        let out = inl - inr;
        out.value().into()
    };
    assert_eq!(exec_out, reference_out);
});

test_binop!(mul, W, |inl, inr, m| {
    let exec_out: W = {
        let executor = ExecutionBackend::<WP>::new().into_executor();
        let inl = MontgomeryWordRef::new(executor.input(inl), m);
        let inr = MontgomeryWordRef::new(executor.input(inr), m);
        let out = inl * inr;
        executor.output(out.value());
        executor.finalize().as_vec()[0]
    };
    let reference_out = {
        let inl = MontgomeryWord::new(inl, m);
        let inr = MontgomeryWord::new(inr, m);
        let out = inl * inr;
        out.value().into()
    };
    assert_eq!(exec_out, reference_out);
});
