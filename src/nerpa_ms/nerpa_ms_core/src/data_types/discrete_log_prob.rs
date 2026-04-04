use crate::data_types::common_types::LogProb;

const MAX_LOG_PROB: LogProb = 0.0; // probability 1
const MIN_LOG_PROB: LogProb = -50.0; // values below are treated as effectively zero probability

/// Largest discrete value we can represent (inclusive range is `0..=MAX_DISCRETE_LOG_PROB`).
const MAX_DISCRETE_LOG_PROB: usize = (1 << 14) - 1; // 16383 => 16384 bits => exactly 256 u64 words
const SCALING_FACTOR: f64 = (MAX_DISCRETE_LOG_PROB as f64) / (MAX_LOG_PROB - MIN_LOG_PROB);

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DiscreteLogProb(pub usize);

impl DiscreteLogProb {
    pub fn logprob_to_centered_discrete_lp(lp: LogProb) -> i32 {
        (lp * SCALING_FACTOR).round() as i32
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
/// - shift_left(lp): add `lp` to all represented values (dropping out-of-range)
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
	    lps.into_iter()
	    .map(DiscreteLogProb::from_logprob)
	    .collect();
	Self::from_dlp_vec(dlps)
    }

    /// Adds `lp` to all discrete log-prob values represented in this set.
    ///
    /// This is implemented as a bit shift by `round(lp * SCALING_FACTOR)`.
    /// Negative shifts move bits to *lower* indices (dropping anything < 0).
    pub fn add_to_all(&self, lp: LogProb) -> DiscreteLogProbSet {
        let delta: i32 = DiscreteLogProb::logprob_to_centered_discrete_lp(lp);
        assert!(delta < 0, "add_to_all: positive log prob not allowed");
        self.shift_towards_zero(delta.unsigned_abs() as usize)
    }

    fn shift_towards_zero(&self, delta: usize) -> DiscreteLogProbSet {
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

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: LogProb, b: LogProb, eps: LogProb) -> bool {
        (a - b).abs() <= eps
    }

    #[test]
    fn dlp_from_logprob_clamps_endpoints() {
        assert_eq!(DiscreteLogProb::from_logprob(MIN_LOG_PROB - 1.0), DiscreteLogProb(0));
        assert_eq!(DiscreteLogProb::from_logprob(MIN_LOG_PROB), DiscreteLogProb(0));

        assert_eq!(
            DiscreteLogProb::from_logprob(MAX_LOG_PROB),
            DiscreteLogProb(MAX_DISCRETE_LOG_PROB)
        );
        assert_eq!(
            DiscreteLogProb::from_logprob(MAX_LOG_PROB + 1.0),
            DiscreteLogProb(MAX_DISCRETE_LOG_PROB)
        );
    }

    #[test]
    fn dlp_roundtrip_is_reasonable() {
        // Not exact because of rounding, but should be within about half a bin.
        let half_bin = 0.5 / SCALING_FACTOR;

        let samples = [
            MIN_LOG_PROB,
            -49.9,
            -25.0,
            -12.345,
            -1.0,
            -0.001,
            MAX_LOG_PROB,
        ];

        for &lp in &samples {
            let d = DiscreteLogProb::from_logprob(lp);
            let lp2 = d.to_logprob();

            // For endpoints, to_logprob returns MIN_LOG_PROB + d/scale, so it matches exactly.
            assert!(
                approx_eq(lp, lp2, half_bin + 1e-12),
                "lp={lp} d={d:?} lp2={lp2} half_bin={half_bin}"
            );

            // Always in range.
            assert!(d.0 <= MAX_DISCRETE_LOG_PROB);
        }
    }

    #[test]
    fn set_union_works() {
        let a = DiscreteLogProbSet::from_dlp_vec(vec![
            DiscreteLogProb(0),
            DiscreteLogProb(10),
            DiscreteLogProb(1000),
        ]);
        let b = DiscreteLogProbSet::from_dlp_vec(vec![
            DiscreteLogProb(10),
            DiscreteLogProb(11),
            DiscreteLogProb(MAX_DISCRETE_LOG_PROB),
        ]);

        let u = a.union(&b);

        // We don't have a `contains` method, so check via iter_desc output.
        let got: Vec<usize> = u.iter_desc().map(|d| d.0).collect();

        // Should contain 0,10,11,1000,MAX.
        assert!(got.contains(&0));
        assert!(got.contains(&10));
        assert!(got.contains(&11));
        assert!(got.contains(&1000));
        assert!(got.contains(&MAX_DISCRETE_LOG_PROB));
        assert_eq!(got.len(), 5);
    }

    #[test]
    fn shift_towards_zero_bitshift_0_is_safe_and_correct() {
        // Use values that make delta a multiple of 64 (bit_shift == 0).
        let s = DiscreteLogProbSet::from_dlp_vec(vec![DiscreteLogProb(0), DiscreteLogProb(64), DiscreteLogProb(128)]);
        let shifted = s.shift_towards_zero(64);

        // Expected: 64->0, 128->64, 0 falls off
        let got: Vec<usize> = shifted.iter_desc().map(|d| d.0).collect();
        assert_eq!(got, vec![64, 0]);
    }

    #[test]
    fn shift_towards_zero_nontrivial_bitshift() {
        // Delta = 1 crosses within a word; easy to sanity check.
        let s = DiscreteLogProbSet::from_dlp_vec(vec![DiscreteLogProb(1), DiscreteLogProb(63), DiscreteLogProb(64)]);
        let shifted = s.shift_towards_zero(1);

        // Expected: 1->0, 63->62, 64->63
        let got: Vec<usize> = shifted.iter_desc().map(|d| d.0).collect();
        assert_eq!(got, vec![63, 62, 0]);
    }

    #[test]
    fn shift_towards_zero_large_delta_erases_all() {
        let s = DiscreteLogProbSet::from_dlp_vec(vec![DiscreteLogProb(0), DiscreteLogProb(123), DiscreteLogProb(MAX_DISCRETE_LOG_PROB)]);
        let shifted = s.shift_towards_zero(MAX_DISCRETE_LOG_PROB + 1);
        assert_eq!(shifted.iter_desc().next(), None);
    }

    #[test]
    fn iter_desc_is_descending_and_unique() {
        let s = DiscreteLogProbSet::from_dlp_vec(vec![
            DiscreteLogProb(5),
            DiscreteLogProb(5), // duplicate on purpose
            DiscreteLogProb(6),
            DiscreteLogProb(100),
            DiscreteLogProb(64),
        ]);

        let got: Vec<usize> = s.iter_desc().map(|d| d.0).collect();

        // Descending order.
        assert_eq!(got, vec![100, 64, 6, 5]);

        // Also implicitly checks no duplicates in iteration result.
    }

    #[test]
    fn add_to_all_panics_on_non_negative_delta() {
        // lp = 0 => delta = 0, which violates assert!(delta < 0)
        let s = DiscreteLogProbSet::from_dlp_vec(vec![DiscreteLogProb(123)]);

        let result = std::panic::catch_unwind(|| {
            let _ = s.add_to_all(0.0);
        });

        assert!(result.is_err());
    }

    #[test]
    fn from_logprob_vec_matches_from_dlp_vec() {
        // Pick a few logprobs and compare to explicit discretization.
        let lps = vec![-50.0, -25.0, -1.0, 0.0];
        let a = DiscreteLogProbSet::from_logprob_vec(lps.clone());

        let dlps: Vec<DiscreteLogProb> = lps.into_iter().map(DiscreteLogProb::from_logprob).collect();
        let b = DiscreteLogProbSet::from_dlp_vec(dlps);

        assert_eq!(a, b);
    }
}
