use nerpa_ms_core::data_types::common_types::LogProb;
use nerpa_ms_core::data_types::discrete_log_prob::DiscreteLogProbSet;
use nerpa_ms_core::data_types::dp_table::DP_Table;

#[test]
fn union_cell_into_cell_is_dst_or_equals_src() {
    let mut dp = DP_Table::new(2, 0, 2);

    let src = dp.idx(0, 0, 0);
    let dst = dp.idx(1, 0, 1);

    *dp.get_mut(0, 0, 0) = DiscreteLogProbSet::from_logprob_vec(vec![0.0]);
    assert!(dp.get(1, 0, 1).is_empty());

    dp.union_cell_into_cell(src, dst);

    assert_eq!(
        dp.get(1, 0, 1),
        &DiscreteLogProbSet::from_logprob_vec(vec![0.0])
    );
}
