from pathlib import Path
import polars as pl
import yaml

nerpa_dir = Path(__file__).parent.parent

mibig_info_path = nerpa_dir / 'data/for_training_and_testing/mibig_bgcs_info.tsv'
mibig_info = pl.read_csv(mibig_info_path, separator='\t')

# mibig_info header: bgc_id  num_a_domains   origins bigscape_families       compound_id     num_recognized_nodes    iso_class_idx    in_approved_matches     deprecated in mibig4

pnrpdb2_info_path = nerpa_dir / 'data/for_training_and_testing/pnrpdb2_info.tsv'
pnrpdb2_info = pl.read_csv(pnrpdb2_info_path, separator='\t')
# pnrpdb2_info header: compound_id     total_num_monomers      num_rban_recognized_monomers    num_nrps_monomers       nrp_variant_iso_class_representative     rban_graph_iso_class_representative


# q: print total number of unique bgcs
unique_bgcs = mibig_info.select(pl.col('bgc_id')).unique()
print(f"Total number of unique BGCs: {unique_bgcs.height}")

# q: print the number of unique BGCs with num_a_domains > 2 and num_nrps_monomers > 2
# for at least one corresponding compound in pnrpdb2_info
bgcs_filtered = (
    mibig_info.join(
        pnrpdb2_info.select(['compound_id', 'num_nrps_monomers']),
        on='compound_id',
        how='inner',
    )
    .filter((pl.col('num_a_domains') > 2) & (pl.col('num_nrps_monomers') > 2))
    .select('bgc_id')
    .unique()
)

print(
    f"Number of unique BGCs with num_a_domains > 2 and num_nrps_monomers > 2: {bgcs_filtered.height}"
)

# q: print the number of unique BGCs with in_approved_matches = True
bgcs_approved_matches = mibig_info.filter(pl.col('in_approved_matches') == True).select('bgc_id').unique()
print(f"Number of unique BGCs with in_approved_matches = True: {bgcs_approved_matches.height}")

# q: print the number of BGCs from bgcs_filtered that don't have bigscape_families appearing in any BGC with in_approved_matches = True. bigscape_families is a simple string, not a list, and cannot be null. We consider a bigscape family to be the same if the strings are exactly the same (case-sensitive).
approved_bigscape_families = (
    mibig_info.filter(pl.col('in_approved_matches') == True)
    .select('bigscape_families')
    .unique()
    .get_column('bigscape_families')
)

bgcs_filtered_without_approved_families = (
    bgcs_filtered.join(
        mibig_info.select(['bgc_id', 'bigscape_families']).unique(),
        on='bgc_id',
        how='left',
    )
    .filter(~pl.col('bigscape_families').is_in(approved_bigscape_families))
    .select('bgc_id')
    .unique()
)

print(
    "Number of BGCs from bgcs_filtered whose bigscape_families do not appear in any BGC with in_approved_matches = True: "
    f"{bgcs_filtered_without_approved_families.height}"
)

# q: print the total number of compounds with num_rban_recognized_monomers > 2 associated with the BGCs from bgcs_filtered_without_approved_families
compounds_rban_gt2 = (
    bgcs_filtered_without_approved_families.join(
        mibig_info.select(['bgc_id', 'compound_id']).unique(),
        on='bgc_id',
        how='inner',
    )
    .join(
        pnrpdb2_info.select(['compound_id', 'num_rban_recognized_monomers']),
        on='compound_id',
        how='inner',
    )
    .filter(pl.col('num_rban_recognized_monomers') > 2)
    .select('compound_id')
    .unique()
)

print(
    "Total number of compounds with num_rban_recognized_monomers > 2 associated with the BGCs from bgcs_filtered_without_approved_families: "
    f"{compounds_rban_gt2.height}"
)
