use nerpa_ms_core::data_types::discrete_log_prob::{
    DiscreteLogProb,
    DiscreteLogProbSet,
    MIN_LOG_PROB,
    MAX_LOG_PROB,
    MAX_DISCRETE_LOG_PROB,
    SCALING_FACTOR
};
use nerpa_ms_core::data_types::common_types::LogProb;
use std::collections::HashSet;

fn approx_eq(a: LogProb, b: LogProb, eps: LogProb) -> bool {
    (a - b).abs() <= eps
}

#[test]
fn dlp_from_logprob_clamps_endpoints() {
    assert_eq!(
        DiscreteLogProb::from_logprob(MIN_LOG_PROB - 1.0),
        DiscreteLogProb(0)
    );
    assert_eq!(
        DiscreteLogProb::from_logprob(MIN_LOG_PROB),
        DiscreteLogProb(0)
    );

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
    let s = DiscreteLogProbSet::from_dlp_vec(vec![
        DiscreteLogProb(0),
        DiscreteLogProb(64),
        DiscreteLogProb(128),
    ]);
    let shifted = s.shift_towards_zero(64);

    // Expected: 64->0, 128->64, 0 falls off
    let got: Vec<usize> = shifted.iter_desc().map(|d| d.0).collect();
    assert_eq!(got, vec![64, 0]);
}

#[test]
fn shift_towards_zero_nontrivial_bitshift() {
    // Delta = 1 crosses within a word; easy to sanity check.
    let s = DiscreteLogProbSet::from_dlp_vec(vec![
        DiscreteLogProb(1),
        DiscreteLogProb(63),
        DiscreteLogProb(64),
    ]);
    let shifted = s.shift_towards_zero(1);

    // Expected: 1->0, 63->62, 64->63
    let got: Vec<usize> = shifted.iter_desc().map(|d| d.0).collect();
    assert_eq!(got, vec![63, 62, 0]);
}

#[test]
fn shift_towards_zero_large_delta_erases_all() {
    let s = DiscreteLogProbSet::from_dlp_vec(vec![
        DiscreteLogProb(0),
        DiscreteLogProb(123),
        DiscreteLogProb(MAX_DISCRETE_LOG_PROB),
    ]);
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
fn add_to_all_panics_on_positive_delta() {
    // lp = 0.1 => delta = 0, which violates assert!(delta < 0)
    let s = DiscreteLogProbSet::from_dlp_vec(vec![DiscreteLogProb(123)]);

    let result = std::panic::catch_unwind(|| {
        let _ = s.add_to_all(0.1);
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
#[test]
fn shift_towards_zero_matches_naive_hashset_many_values_many_shifts() {
    // Build a fairly dense set: every 5th value, plus some "interesting" boundaries.
    let mut src_vals: HashSet<usize> = (0..=MAX_DISCRETE_LOG_PROB).step_by(5).collect();
    src_vals.extend([1, 2, 3, 4, 63, 64, 65, 127, 128, 129, MAX_DISCRETE_LOG_PROB]);

    let s =
        DiscreteLogProbSet::from_dlp_vec(src_vals.iter().copied().map(DiscreteLogProb).collect());

    // Try lots of deltas, including word-boundary-ish ones and very large ones.
    let mut deltas: Vec<usize> = (0..=300).collect();
    deltas.extend([
        63, 64, 65, 127, 128, 129, 255, 256, 257, 511, 512, 513, 1000, 4096,
    ]);
    deltas.extend([
        MAX_DISCRETE_LOG_PROB - 1,
        MAX_DISCRETE_LOG_PROB,
        MAX_DISCRETE_LOG_PROB + 1,
        MAX_DISCRETE_LOG_PROB + 1000,
    ]);

    for delta in deltas {
        // Naive expected behavior:
        // shifting "towards zero" by delta maps v -> v - delta, dropping anything < 0.
        let mut expected: HashSet<usize> = HashSet::new();
        for &v in &src_vals {
            if v >= delta {
                expected.insert(v - delta);
            }
        }

        let got: HashSet<usize> = s
            .shift_towards_zero(delta)
            .iter_desc()
            .map(|d| d.0)
            .collect();

        assert_eq!(got, expected, "mismatch for delta={delta}");
    }
}
