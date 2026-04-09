use std::collections::HashSet;

use crate::data_types::common_types::LogProb;

pub const MAX_LOG_PROB: LogProb = 0.0; // probability 1
pub const MIN_LOG_PROB: LogProb = -50.0; // values below are treated as effectively zero probability

/// Largest discrete value we can represent (inclusive range is `0..=MAX_DISCRETE_LOG_PROB`).
pub const MAX_DISCRETE_LOG_PROB: usize = (1 << 14) - 1; // 16383 => 16384 bits => exactly 256 u64 words
pub const SCALING_FACTOR: f64 = (MAX_DISCRETE_LOG_PROB as f64) / (MAX_LOG_PROB - MIN_LOG_PROB);

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DiscreteLogProb(pub usize);

impl DiscreteLogProb {
    pub fn logprob_to_centered_discrete_lp(lp: LogProb) -> i64 {
        (lp * SCALING_FACTOR).round() as i64
    }

    pub fn from_logprob(lp: LogProb) -> Self {
        if lp <= MIN_LOG_PROB {
            DiscreteLogProb(0)
        } else if lp >= MAX_LOG_PROB {
            DiscreteLogProb(MAX_DISCRETE_LOG_PROB)
        } else {
            let d = ((lp - MIN_LOG_PROB) * SCALING_FACTOR).round() as usize;
            DiscreteLogProb(d)
        }
    }

    pub fn to_logprob(self) -> LogProb {
        MIN_LOG_PROB + (self.0 as f64) / SCALING_FACTOR
    }

    pub fn shift(self, delta: i64) -> DiscreteLogProb {
	let new_d = (self.0 as i64) + delta;
	if new_d < 0 || new_d > MAX_DISCRETE_LOG_PROB as i64 {
	    panic!("shift: resulting discrete log prob out of range: new_d={new_d}");
	} else {
	    DiscreteLogProb(new_d as usize)
	}
    }
}

impl From<DiscreteLogProb> for usize {
    fn from(dlp: DiscreteLogProb) -> Self {
        dlp.0
    }
}

/// Number of 64-bit words needed to store `MAX_DISCRETE_LOG_PROB + 1` bits.
/// With MAX = 2^k - 1 and k >= 6, this is an integer and divisible by 64.
const N_WORDS: usize = (MAX_DISCRETE_LOG_PROB + 1) / 64;

// Compile-time sanity check: number of bits must be divisible by 64.
const _: [(); 0] = [(); (MAX_DISCRETE_LOG_PROB + 1) % 64];

/// Bit-array-like set of discrete log-prob values in `0..=MAX_DISCRETE_LOG_PROB`.
///
/// Only supports what you asked for:
/// - shift_towards_zero(lp): add `lp` to all represented values (dropping out-of-range)
/// - union(other): bitwise OR
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DiscreteLogProbSet {
    words: [u64; N_WORDS],
}

impl DiscreteLogProbSet {
    pub fn empty() -> Self {
        Self {
            words: [0u64; N_WORDS],
        }
    }

    pub fn is_empty(&self) -> bool {
	self.words.iter().all(|&w| w == 0)
    }

    pub fn contains(&self, dlp: DiscreteLogProb) -> bool {
	let i: usize = dlp.into();
	(self.words[i / 64] & (1u64 << (i % 64))) != 0
    }

    pub fn from_dlp_vec(dlps: Vec<DiscreteLogProb>) -> Self {
        let mut s = Self::empty();
        for dlp in dlps {
            let i: usize = dlp.into();
            s.words[i / 64] |= 1u64 << (i % 64);
        }
        s
    }

    pub fn from_logprob_vec(lps: Vec<LogProb>) -> Self {
        let dlps: Vec<DiscreteLogProb> =
            lps.into_iter().map(DiscreteLogProb::from_logprob).collect();
        Self::from_dlp_vec(dlps)
    }

    pub fn get_abs_shift(lp: LogProb) -> usize {
	let delta: i64 = DiscreteLogProb::logprob_to_centered_discrete_lp(lp);
	assert!(delta <= 0, "get_abs_shift: positive log prob not allowed");
	delta.unsigned_abs() as usize
    }

    /// Adds `lp` to all discrete log-prob values represented in this set.
    ///
    /// This is implemented as a bit shift by `round(lp * SCALING_FACTOR)`.
    /// Negative shifts move bits to *lower* indices (dropping anything < 0).
    pub fn add_to_all(&self, lp: LogProb) -> DiscreteLogProbSet {
        let delta = DiscreteLogProbSet::get_abs_shift(lp);
        self.shift_towards_zero(delta)
    }

    pub fn shift_towards_zero(&self, delta: usize) -> DiscreteLogProbSet {
        if delta == 0 {
            return self.clone();
        }

        let mut out = DiscreteLogProbSet::empty();

        let (word_shift, bit_shift) = (delta / 64, delta % 64);

        for dst in 0..N_WORDS.saturating_sub(word_shift) {
            // preimage of bits [dst*64 .. (dst+1)*64] is bits [dst*64 + delta .. (dst+1)*64 + delta]:
            // `dst + word_shift` bits [bit_shift .. 63] map to src[0: 64 - bit_shift], and
            // `dst + word_shift + 1` bits [0 .. bit_shift-1] map to src[64 - bit_shift : 64]
            let src1 = self.words[dst + word_shift];
            let shifted_src1 = src1 >> bit_shift;

            let src2: u64 = if dst + word_shift + 1 < N_WORDS {
                self.words[dst + word_shift + 1]
            } else {
                0
            };

            // shifting by 64 or more is UB in Rust
            let shifted_src2 = if bit_shift > 0 {
                src2 << (64 - bit_shift)
            } else {
                0
            };

            out.words[dst] = shifted_src1 | shifted_src2;
        }

        out
    }

    /// Bitwise OR union.
    pub fn union(&self, other: &DiscreteLogProbSet) -> DiscreteLogProbSet {
        let mut out = self.clone();
        for i in 0..N_WORDS {
            out.words[i] |= other.words[i];
        }
        out
    }

    /// Bitwise OR union.
    pub fn union_inplace(&mut self, other: &DiscreteLogProbSet) -> () {
        for i in 0..N_WORDS {
            self.words[i] |= other.words[i];
        }
    }

    /// Iterate over contained values from largest to smallest.
    pub fn iter_desc(&self) -> DiscreteLogProbSetIterDesc<'_> {
        let mut it = DiscreteLogProbSetIterDesc {
            set: self,
            word_idx: N_WORDS as isize - 1,
            cur_word: 0,
        };

        if it.word_idx >= 0 {
            it.cur_word = self.words[it.word_idx as usize];
        }
        it
    }
}

pub struct DiscreteLogProbSetIterDesc<'a> {
    set: &'a DiscreteLogProbSet,
    word_idx: isize,
    cur_word: u64,
}

impl<'a> Iterator for DiscreteLogProbSetIterDesc<'a> {
    type Item = DiscreteLogProb;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.word_idx < 0 {
                return None;
            }

            if self.cur_word != 0 {
                // Take the highest set bit in the current word.
                let bit = 63usize - (self.cur_word.leading_zeros() as usize);
                self.cur_word &= !(1u64 << bit);

                let idx: usize = (self.word_idx as usize) * 64 + bit;
                return Some(DiscreteLogProb(idx));
            }

            // Move to the next lower word.
            self.word_idx -= 1;
            if self.word_idx >= 0 {
                self.cur_word = self.set.words[self.word_idx as usize];
            }
        }
    }
}
