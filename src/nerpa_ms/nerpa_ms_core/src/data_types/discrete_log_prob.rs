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

    pub fn singleton_from_dlp(dlp: DiscreteLogProb) -> Self {
        let mut s = Self::empty();
        let i: usize = dlp.into();
        if i <= MAX_DISCRETE_LOG_PROB {
            s.words[i / 64] |= 1u64 << (i % 64);
        }
        s
    }

    pub fn singleton_from_lp(lp: LogProb) -> Self {
	Self::singleton_from_dlp(DiscreteLogProb::from_logprob(lp))
    }

    /// Adds `lp` to all discrete log-prob values represented in this set.
    ///
    /// This is implemented as a bit shift by `round(lp * SCALING_FACTOR)`.
    /// Negative shifts move bits to *lower* indices (dropping anything < 0).
    pub fn add_to_all(&self, lp: LogProb) -> DiscreteLogProbSet {
        let delta: i32 = DiscreteLogProb::logprob_to_centered_discrete_lp(lp);
        assert!(delta < 0, "add_to_all: positive log prob not allowed");
        self.shift_left(delta.unsigned_abs() as usize)
    }

    fn shift_left(&self, delta: usize) -> DiscreteLogProbSet {
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
            let shifted_src1 = src1 << bit_shift;

            let src2: u64 = if dst + word_shift + 1 < N_WORDS {
                self.words[dst + word_shift + 1]
            } else {
                0
            };

            // shifting by 64 or more is UB in Rust
            let shifted_src2 = if bit_shift > 0 {
                src2 >> (64 - bit_shift)
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
}
