// SPDX-License-Identifier: LGPL-3.0-or-later

//! Head-to-head gate-count comparison of the two secp256k1 field-reduction strategies:
//! Montgomery (with NAF constant-multiply) vs pseudo-Mersenne (Solinas) reduction.
//!
//! Run with the local zkboo tree:
//! ```text
//! cargo run --release --example profile_reduction \
//!   --config 'patch.crates-io.zkboo.path="../zkboo"' \
//!   --config 'patch.crates-io.zkboo-profiling.path="../zkboo-profiling"'
//! ```

use zkboo::backend::{Backend, Frontend};
use zkboo::circuit::Circuit;
use zkboo::word::CompositeWord;
use zkboo_modular::montgomery::{MontgomeryMod, MontgomeryWordRef};
use zkboo_modular::pseudo_mersenne::{self, PseudoMersenneMod};
use zkboo_profiling::profile;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SecpMont;
impl MontgomeryMod<u64, 4> for SecpMont {
    fn n(&self) -> CompositeWord<u64, 4> {
        CompositeWord::from_be_words([
            0xffffffffffffffff,
            0xffffffffffffffff,
            0xffffffffffffffff,
            0xfffffffefffffc2f,
        ])
    }
    fn rr_mod_n(&self) -> CompositeWord<u64, 4> {
        CompositeWord::from_be_words([0, 0, 1, 0x000007a2000e90a1])
    }
    fn n_neg_inv(&self) -> CompositeWord<u64, 4> {
        CompositeWord::from_be_words([
            0xc9bd190515538399,
            0x9c46c2c295f2b761,
            0xbcb223fedc24a059,
            0xd838091dd2253531,
        ])
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SecpPM;
impl PseudoMersenneMod<u64, 4> for SecpPM {
    fn p(&self) -> CompositeWord<u64, 4> {
        SecpMont.n()
    }
    fn c(&self) -> CompositeWord<u64, 4> {
        CompositeWord::from_be_words([0, 0, 0, 0x0000000100_0003d1])
    }
}

fn w<B: Backend>(fe: &Frontend<B>) -> zkboo::backend::WordRef<B, u64, 4> {
    fe.input(CompositeWord::from_le_words([3, 5, 7, 9]))
}

struct MontMul;
impl Circuit for MontMul {
    fn exec<B: Backend>(&self, fe: &Frontend<B>) {
        let a = MontgomeryWordRef::from_inner(w(fe), SecpMont);
        let b = MontgomeryWordRef::from_inner(w(fe), SecpMont);
        fe.output((a * b).into_inner());
    }
}
/// Montgomery multiply including the input/output domain conversions (to/from Montgomery),
/// which pseudo-Mersenne does not need.
struct MontMulWithConv;
impl Circuit for MontMulWithConv {
    fn exec<B: Backend>(&self, fe: &Frontend<B>) {
        let a = MontgomeryWordRef::new(w(fe), SecpMont);
        let b = MontgomeryWordRef::new(w(fe), SecpMont);
        fe.output((a * b).value());
    }
}
struct PmMul;
impl Circuit for PmMul {
    fn exec<B: Backend>(&self, fe: &Frontend<B>) {
        fe.output(pseudo_mersenne::mul(&SecpPM, w(fe), w(fe)));
    }
}
struct PmReduce;
impl Circuit for PmReduce {
    fn exec<B: Backend>(&self, fe: &Frontend<B>) {
        fe.output(pseudo_mersenne::reduce_wide(&SecpPM, w(fe), w(fe)));
    }
}

fn report<C: Circuit>(label: &str, c: &C) {
    println!("{label:<34} nl_and_msgs = {:>8}", profile(c).and_msg_size().sum());
}

fn main() {
    println!("=== secp256k1 reduction strategy comparison ===");
    report("Montgomery mul (no conversion)", &MontMul);
    report("Montgomery mul (+ to/from conv)", &MontMulWithConv);
    report("pseudo-Mersenne reduce_wide", &PmReduce);
    report("pseudo-Mersenne mul", &PmMul);
}
