use nerpa_ms_core::data_types::discrete_log_prob::{
    DiscreteLogOdds,
    DiscreteLogOddsSet,
    MIN_LOG_PROB,
    MAX_LOG_PROB,
    MAX_DISCRETE_LOG_PROB,
    SCALING_FACTOR
};
use nerpa_ms_core::data_types::common_types::LogOdds;
use std::collections::HashSet;

fn approx_eq(a: LogOdds, b: LogOdds, eps: LogOdds) -> bool {
    (a - b).abs() <= eps
}

#[test]
fn dlp_from_logprob_clamps_endpoints() {
    assert_eq!(
        DiscreteLogOdds::from_logodds(MIN_LOG_PROB - 1.0),
        DiscreteLogOdds(0)
    );
    assert_eq!(
        DiscreteLogOdds::from_logodds(MIN_LOG_PROB),
        DiscreteLogOdds(0)
    );

    assert_eq!(
        DiscreteLogOdds::from_logodds(MAX_LOG_PROB),
        DiscreteLogOdds(MAX_DISCRETE_LOG_PROB)
    );
    let result = std::panic::catch_unwind(|| {
        let _ = DiscreteLogOdds::from_logodds(MAX_LOG_PROB + 1.0);
    });

    assert!(result.is_err());
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
        let d = DiscreteLogOdds::from_logodds(lp);
        let lp2 = d.to_logodds();

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
    let a = DiscreteLogOddsSet::from_dlo_vec(vec![
        DiscreteLogOdds(0),
        DiscreteLogOdds(10),
        DiscreteLogOdds(1000),
    ]);
    let b = DiscreteLogOddsSet::from_dlo_vec(vec![
        DiscreteLogOdds(10),
        DiscreteLogOdds(11),
        DiscreteLogOdds(MAX_DISCRETE_LOG_PROB),
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
    let s = DiscreteLogOddsSet::from_dlo_vec(vec![
        DiscreteLogOdds(0),
        DiscreteLogOdds(64),
        DiscreteLogOdds(128),
    ]);
    let shifted = s.shift(-64);

    // Expected: 64->0, 128->64, 0 falls off
    let got: Vec<usize> = shifted.iter_desc().map(|d| d.0).collect();
    assert_eq!(got, vec![64, 0]);
}

#[test]
fn shift_towards_zero_nontrivial_bitshift() {
    // Delta = 1 crosses within a word; easy to sanity check.
    let s = DiscreteLogOddsSet::from_dlo_vec(vec![
        DiscreteLogOdds(1),
        DiscreteLogOdds(63),
        DiscreteLogOdds(64),
    ]);
    let shifted = s.shift(1);

    // Expected: 1->0, 63->62, 64->63
    let got: Vec<usize> = shifted.iter_desc().map(|d| d.0).collect();
    assert_eq!(got, vec![63, 62, 0]);
}

#[test]
fn shift_towards_zero_large_delta_erases_all() {
    let s = DiscreteLogOddsSet::from_dlo_vec(vec![
        DiscreteLogOdds(0),
        DiscreteLogOdds(123),
        DiscreteLogOdds(MAX_DISCRETE_LOG_PROB),
    ]);
    let shifted = s.shift(-(MAX_DISCRETE_LOG_PROB as isize) - 1);
    assert_eq!(shifted.iter_desc().next(), None);
}

#[test]
fn iter_desc_is_descending_and_unique() {
    let s = DiscreteLogOddsSet::from_dlo_vec(vec![
        DiscreteLogOdds(5),
        DiscreteLogOdds(5), // duplicate on purpose
        DiscreteLogOdds(6),
        DiscreteLogOdds(100),
        DiscreteLogOdds(64),
    ]);

    let got: Vec<usize> = s.iter_desc().map(|d| d.0).collect();

    // Descending order.
    assert_eq!(got, vec![100, 64, 6, 5]);

    // Also implicitly checks no duplicates in iteration result.
}

#[test]
fn add_to_all_panics_on_positive_overflow() {
    let s = DiscreteLogOddsSet::from_dlo_vec(vec![DiscreteLogOdds(MAX_DISCRETE_LOG_PROB)]);

    let result = std::panic::catch_unwind(|| {
        let _ = s.add_to_all(1.0);
    });

    assert!(result.is_err());
}

#[test]
fn from_logprob_vec_matches_from_dlp_vec() {
    // Pick a few logprobs and compare to explicit discretization.
    let lps = vec![-50.0, -25.0, -1.0, 0.0];
    let a = DiscreteLogOddsSet::from_logodds_vec(lps.clone());

    let dlps: Vec<DiscreteLogOdds> = lps.into_iter().map(DiscreteLogOdds::from_logodds).collect();
    let b = DiscreteLogOddsSet::from_dlo_vec(dlps);

    assert_eq!(a, b);
}
#[test]
fn shift_matches_naive_hashset_many_values_many_shifts() {
    // Build a fairly dense set: every 5th value, plus some "interesting" boundaries.
    let mut src_vals: HashSet<usize> = (0..=MAX_DISCRETE_LOG_PROB).step_by(5).collect();
    src_vals.extend([1, 2, 3, 4, 63, 64, 65, 127, 128, 129, MAX_DISCRETE_LOG_PROB]);

    let s =
        DiscreteLogOddsSet::from_dlo_vec(src_vals.iter().copied().map(DiscreteLogOdds).collect());

    // Try lots of deltas, including word-boundary-ish ones and very large ones.
    let deltas_pos: Vec<isize> = {
	let mut deltas_pos: Vec<isize> = (0..=300).collect();
	deltas_pos.extend([
            63, 64, 65, 127, 128, 129, 255, 256, 257, 511, 512, 513, 1000, 4096,
	]);
	deltas_pos.extend([
            (MAX_DISCRETE_LOG_PROB as isize) - 1,
            (MAX_DISCRETE_LOG_PROB as isize),
            (MAX_DISCRETE_LOG_PROB as isize) + 1,
            (MAX_DISCRETE_LOG_PROB as isize) + 1000,
	]);
	deltas_pos
    };
    let deltas_neg: Vec<isize> = {
	deltas_pos.iter().map(|&d| -d).collect()
    };

    for delta in deltas_neg {
        // Naive expected behavior:
        // shifting "towards zero" by delta maps v -> v - delta, dropping anything < 0.
        let mut expected: Vec<usize> = Vec::new();
        for &v in &src_vals {
            if (v as isize) + delta >= 0 {
                expected.push(((v as isize) + delta) as usize);
            }
        }
	expected.sort_unstable();

        let mut got: Vec<usize> = s
            .shift(delta)
            .iter_desc()
            .map(|d| d.0)
            .collect();

	got.sort_unstable();

	let first_mismatch: Option<(Option<usize>, Option<usize>)> = {
	    let mut res = None;
	    for idx in 0..got.len() {
		if idx >= expected.len() {
		    res = Some((None, Some(got[idx])));
		    break;
		} else if got[idx] != expected[idx] {
		    res = Some((Some(expected[idx]), Some(got[idx])));
		    break;
		}
	    }
	    if res.is_none() && expected.len() > got.len() {
		res = Some((Some(expected[got.len()]), None));
	    }
	    res
	};

	assert!(first_mismatch.is_none(), "mismatch for delta={delta}\nFirst mismatch at expected={:?}, got={:?}", first_mismatch.unwrap().0, first_mismatch.unwrap().1);
    }

    for delta in deltas_pos {
	// positive overflow panics, so skip any deltas that would cause overflow for any value in src_vals.
	let src_vals_filtered: Vec<usize> = {
	    src_vals.iter()
		.filter(|&&v| (v as isize) + delta <= MAX_DISCRETE_LOG_PROB as isize)
		.copied()
		.collect()
	};
	let s = {
	    let dlo_vec: Vec<DiscreteLogOdds> = {
		src_vals_filtered.iter()
		    .copied()
		    .map(DiscreteLogOdds)
		    .collect()
	    };
            DiscreteLogOddsSet::from_dlo_vec(dlo_vec)
	};

	let mut expected: Vec<usize> = Vec::new();
	for &v in &src_vals {
	    expected.push(((v as isize) + delta) as usize);
	}
	expected.sort_unstable();

	let mut got: Vec<usize> = s
	    .shift(delta)
	    .iter_desc()
	    .map(|d| d.0)
	    .collect();

	got.sort_unstable();

	let first_mismatch: Option<(Option<usize>, Option<usize>)> = {
	    let mut res = None;
	    for idx in 0..got.len() {
		if idx >= expected.len() {
		    res = Some((None, Some(got[idx])));
		    break;
		} else if got[idx] != expected[idx] {
		    res = Some((Some(expected[idx]), Some(got[idx])));
		    break;
		}
	    }
	    if res.is_none() && expected.len() > got.len() {
		res = Some((Some(expected[got.len()]), None));
	    }
	    res
	};

	assert!(first_mismatch.is_none(), "mismatch for delta={delta}\nFirst mismatch at expected={:?}, got={:?}", first_mismatch.unwrap().0, first_mismatch.unwrap().1);
    }
}
