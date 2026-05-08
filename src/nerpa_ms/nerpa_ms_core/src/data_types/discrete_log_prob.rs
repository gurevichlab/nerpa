use std::collections::HashSet;

use crate::data_types::common_types::LogOdds;

pub const MAX_LOG_PROB: LogOdds = 100.0; // assumed to be unreachable. Values above this cause overflow and panic
pub const MIN_LOG_PROB: LogOdds = -50.0; // values below are discarded as impossible

/// Largest discrete value we can represent (inclusive range is `0..=MAX_DISCRETE_LOG_PROB`).
pub const MAX_DISCRETE_LOG_PROB: usize = (1 << 14) - 1; // 16383 => 16384 bits => exactly 256 u64 words
pub const SCALING_FACTOR: f64 = (MAX_DISCRETE_LOG_PROB as f64) / (MAX_LOG_PROB - MIN_LOG_PROB);

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DiscreteLogOdds(pub usize);

impl DiscreteLogOdds {
    pub fn logodds_to_centered_discrete_lo(lp: LogOdds) -> isize {
        (lp * SCALING_FACTOR).round() as isize
    }

    pub fn from_logodds(lo: LogOdds) -> Self {
        if lo <= MIN_LOG_PROB {
            DiscreteLogOdds(0)
        } else if lo > MAX_LOG_PROB {
	    panic!("DiscreteLogOdds::from_logodds: log odds {} above MAX_LOG_PROB {}", lo, MAX_LOG_PROB);
        } else {
            let d = ((lo - MIN_LOG_PROB) * SCALING_FACTOR).round() as usize;
            DiscreteLogOdds(d)
        }
    }

    pub fn to_logodds(self) -> LogOdds {
        MIN_LOG_PROB + (self.0 as f64) / SCALING_FACTOR
    }

    pub fn shift(self, delta: isize) -> Option<DiscreteLogOdds> {
	let new_d = (self.0 as isize) + delta;
	if new_d < 0 || new_d > MAX_DISCRETE_LOG_PROB as isize {
	    None
	} else {
	    Some(DiscreteLogOdds(new_d as usize))
	}
    }
}

impl From<DiscreteLogOdds> for usize {
    fn from(dlo: DiscreteLogOdds) -> Self {
        dlo.0
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
pub struct DiscreteLogOddsSet {
    words: [u64; N_WORDS],
}

impl DiscreteLogOddsSet {
    pub fn empty() -> Self {
        Self {
            words: [0u64; N_WORDS],
        }
    }

    pub fn is_empty(&self) -> bool {
	self.words.iter().all(|&w| w == 0)
    }

    pub fn contains(&self, dlo: DiscreteLogOdds) -> bool {
	let i: usize = dlo.into();
	(self.words[i / 64] & (1u64 << (i % 64))) != 0
    }

    pub fn from_dlo_vec(dlos: Vec<DiscreteLogOdds>) -> Self {
        let mut s = Self::empty();
        for dlo in dlos {
            let i: usize = dlo.into();
            s.words[i / 64] |= 1u64 << (i % 64);
        }
        s
    }

    pub fn from_logodds_vec(los: Vec<LogOdds>) -> Self {
        let dlps: Vec<DiscreteLogOdds> =
            los.into_iter().map(DiscreteLogOdds::from_logodds).collect();
        Self::from_dlo_vec(dlps)
    }

    pub fn get_abs_shift(lp: LogOdds) -> usize {
	let delta: isize = DiscreteLogOdds::logodds_to_centered_discrete_lo(lp);
	assert!(delta <= 0, "get_abs_shift: positive log prob not allowed");
	delta.unsigned_abs() as usize
    }

    fn check_shift_overflow(&self, delta: isize) -> bool {
	let max_element: Option<usize> = {
	    self.words.iter()
		.rposition(|&w| w != 0)
		.map(|i| i * 64 + (63 - self.words[i].leading_zeros() as usize))
	};
	if let Some(m) = max_element {
	    if delta > 0 && (m as isize) + delta > MAX_DISCRETE_LOG_PROB as isize {
		return true;
	    }
	}

	let min_element: Option<usize> = {
	    self.words.iter()
		.position(|&w| w != 0)
		.map(|i| i * 64 + self.words[i].trailing_zeros() as usize)
	};
	if let Some(m) = min_element {
	    if delta < 0 && (m as isize) + delta < 0 {
		return true;
	    }
	}

	false
    }

    /// Adds `lo` to all discrete log-odds values represented in this set.
    ///
    /// This is implemented as a bit shift by `round(lp * SCALING_FACTOR)`.
    /// Negative shifts move bits to *lower* indices (dropping anything < 0).
    pub fn add_to_all(&self, lo: LogOdds) -> DiscreteLogOddsSet {
	let delta: isize = DiscreteLogOdds::logodds_to_centered_discrete_lo(lo);
        self.shift(delta)
    }

    pub fn shift(&self, delta: isize) -> DiscreteLogOddsSet {
        if delta == 0 {
            return self.clone();
        }
	debug_assert!(delta < 0 || !self.check_shift_overflow(delta),  // overflow for negative shift is fine
		      "shift: overflow detected for delta={}", delta);

        let mut out = DiscreteLogOddsSet::empty();

        let (word_shift, bit_shift) = if delta >= 0 {
	    (delta / 64, delta % 64)
	} else {
	    (delta.div_euclid(64), delta.rem_euclid(64))
	};

	let (dst_lb, dst_ub) = if word_shift >= 0 {  // semi-open: [dst_lb, dst_ub) 
	    (word_shift as usize, N_WORDS)
	} else {
	    (0, N_WORDS.saturating_sub((-word_shift) as usize))
	};

        for dst in dst_lb..dst_ub {
            // preimage of bits [dst*64 .. (dst+1)*64] is bits [dst*64 + delta .. (dst+1)*64 + delta]:
            // `dst + word_shift` bits [bit_shift .. 63] map to src[0: 64 - bit_shift], and
            // `dst + word_shift + 1` bits [0 .. bit_shift-1] map to src[64 - bit_shift : 64]
	    let src1_idx = (dst as isize) - word_shift;
            let src1: u64 = if src1_idx >= 0 && (src1_idx as usize) < N_WORDS {
		self.words[src1_idx as usize]
	    } else {
		0
	    };
            let shifted_src1 = src1 >> bit_shift;

	    let src2_idx = (dst as isize) - word_shift + 1;
            let src2: u64 = if src2_idx >= 0 && (src2_idx as usize) < N_WORDS {
                self.words[src2_idx as usize]
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
    pub fn union(&self, other: &DiscreteLogOddsSet) -> DiscreteLogOddsSet {
        let mut out = self.clone();
        for i in 0..N_WORDS {
            out.words[i] |= other.words[i];
        }
        out
    }

    /// Bitwise OR union.
    pub fn union_inplace(&mut self, other: &DiscreteLogOddsSet) -> () {
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
    set: &'a DiscreteLogOddsSet,
    word_idx: isize,  // signed for convenient loop termination at -1
    cur_word: u64,
}

impl<'a> Iterator for DiscreteLogProbSetIterDesc<'a> {
    type Item = DiscreteLogOdds;

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
                return Some(DiscreteLogOdds(idx));
            }

            // Move to the next lower word.
            self.word_idx -= 1;
            if self.word_idx >= 0 {
                self.cur_word = self.set.words[self.word_idx as usize];
            }
        }
    }
}
