import polars as pl
from pathlib import Path
import argparse
import yaml

from src.antismash_parsing.bgc_variant_types import BGC_Variant


def parse_args() -> argparse.Namespace:
    # q: args -- nerpa_results_dir and output
    parser = argparse.ArgumentParser(description="Add metadata to nerpa report.")
    parser.add_argument(
        "--nerpa-results",
        type=Path,
        required=True,
        help="Directory with nerpa results, expected to contain preprocessed_input/BGC_variants.yaml and the original input TSV with SMILES.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Directory where the enriched report will be written.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    nerpa_dir = Path(__file__).parent.parent
    pnrpdb_path = nerpa_dir / "data/input/pnrpdb2.tsv"

    # 1. Load metadata for BGCs
    bgc_variants_path = args.nerpa_results / 'preprocessed_input' / 'BGC_variants.yaml'
    with open(bgc_variants_path, 'r') as f:
        bgc_variants = [
            BGC_Variant.from_yaml_dict(record)
            for record in yaml.safe_load(f)
        ]

    genux_by_bgc_id = {
        variant.metadata.antismash_json: variant.metadata.taxonomy[-1].strip().lower()
        for variant in bgc_variants
    }

    # 2. Load metadata for NRPs
    pnrpdb = pl.read_csv(pnrpdb_path, separator='\t')
    # q: extract dict "ID" ->"Origin"
    _nrp_origin_by_id = (
        pnrpdb.select(
            pl.col("ID").cast(pl.Utf8).alias("ID"),
            pl.col("Origin")
            .cast(pl.Utf8)
            .str.strip_chars()
            .str.to_lowercase()
            .alias("Origin"),
        )
        .drop_nulls(["ID", "Origin"])
    )
    nrp_origin_by_id = dict(zip(_nrp_origin_by_id["ID"].to_list(),
                                _nrp_origin_by_id["Origin"].to_list()))

    def _origin_to_genera(origin: str | None) -> list[str]:
        # Origin is a ';' separated list of species; we match on genus (first token) for any entry.
        if origin is None:
            return []
        parts: list[str] = [p.strip().lower() for p in origin.split(';') if p.strip()]
        genera: list[str] = [p.split()[0] for p in parts if p]
        return list(dict.fromkeys(genera))  # de-dupe, keep order

    genera_by_nrp_id: dict[str, list[str]] = {
        nrp_id: _origin_to_genera(origin)
        for nrp_id, origin in nrp_origin_by_id.items()
    }

    # q: load report. For each compound, add columns "BGC_genus" and "NRP_genus" using the dicts above. For the IDs, use the report columns "NRP_ID" and "antiSMASH_input"
    report_path = args.nerpa_results / "report.tsv"
    report = pl.read_csv(report_path, separator='\t', infer_schema_length=10000)

    report = report.with_columns(
        pl.col("antiSMASH_input")
        .cast(pl.Utf8)
        .map_elements(lambda x: genux_by_bgc_id.get(x) if x is not None else None, return_dtype=pl.Utf8)
        .alias("BGC_genus"),
        pl.col("NRP_ID")
        .cast(pl.Utf8)
        .map_elements(lambda x: genera_by_nrp_id.get(x) if x is not None else None, return_dtype=pl.List(pl.Utf8))
        .alias("NRP_genera"),
        )

    # q: compute cumulative fraction of rows with BGC_genus == NRP_genus (None counts as not equal). Write this to "cumulative_fraction_genus_match" column. The report is expected to be sorted by score descending, so this will give the cumulative fraction of matches with genus match as we go down the ranked list.
    match_expr = (
         (pl.col("BGC_genus").is_not_null())
        & (pl.col("NRP_genera").is_not_null())
        & (pl.col("NRP_genera").list.contains(pl.col("BGC_genus")))
    ).cast(pl.Int64)
    n_rows: int = report.height
    rank: pl.Series = pl.Series('rank', range(1, n_rows + 1))

 
    report = report.with_columns(
        (match_expr.cum_sum() / rank)
        .cast(pl.Float64)
        .alias("cumulative_fraction_genus_match")
    )

    # CSV/TSV can't store nested list columns; stringify the genera list for output.
    report_for_output = report.with_columns(
        pl.when(pl.col('NRP_genera').is_null())
        .then(None)
        .otherwise(pl.col('NRP_genera').list.join(';'))
        .alias('NRP_genera_str')
    ).drop('NRP_genera').rename({'NRP_genera_str': 'NRP_genera'})

    args.output.mkdir(parents=True, exist_ok=True)
    out_path = args.output / f'{report_path.stem}_with_genus.tsv'
    report_for_output.write_csv(out_path, separator='\t')


    # q: compute probability that a random BGC-NRP pair has the same genus. Compute genus frequencies for each part and then take the scalar product.
    bgc_freq = (
        report.select(pl.col("BGC_genus").cast(pl.Utf8).alias("genus"))
        .drop_nulls(["genus"])
        .group_by("genus")
        .len()
        .with_columns((pl.col("len") / pl.col("len").sum()).alias("bgc_freq"))
        .select(["genus", "bgc_freq"])
    )

    nrp_freq = (
        report.select(pl.col("NRP_genera"))
        .drop_nulls(["NRP_genera"])
        .explode("NRP_genera")
        .drop_nulls(["NRP_genera"])
        .select(pl.col("NRP_genera").cast(pl.Utf8).alias("genus"))
        .group_by("genus")
        .len()
        .with_columns((pl.col("len") / pl.col("len").sum()).alias("nrp_freq"))
        .select(["genus", "nrp_freq"])
    )

    random_pair_same_genus_prob = (
        bgc_freq.join(nrp_freq, on="genus", how="inner")
        .select((pl.col("bgc_freq") * pl.col("nrp_freq")).sum().alias("prob"))
        .item()
    )

    print(f"Random BGC-NRP genus match probability (from marginal genus frequencies): {random_pair_same_genus_prob:.6f}")

    
    # q: plot cumulative_fraction_genus_match vs rank (row number) and save as "cumulative_fraction_genus_match.png" in the output directory
    import matplotlib.pyplot as plt

    n_rows = report.height
    ranks = list(range(1, n_rows + 1))
    y = report["cumulative_fraction_genus_match"].to_list()

    plt.figure(figsize=(8, 5), dpi=150)
    plt.plot(ranks, y, linewidth=2)
    plt.xlabel("Rank (row number)")
    plt.ylabel("Cumulative fraction genus match")
    plt.ylim(0, 1)
    plt.xlim(1, max(1, n_rows))
    plt.grid(True, alpha=0.3)

    # q: add horizontal dashed line for random_pair_same_genus_prob
    plt.axhline(
    y=random_pair_same_genus_prob,
    linestyle="--",
    linewidth=2,
    color="red",
    alpha=0.8,
    label=f"Random pair prob = {random_pair_same_genus_prob:.3f}",
)
    plt.legend(loc="best")

    plot_path = args.output / "cumulative_fraction_genus_match.png"
    plt.tight_layout()
    plt.savefig(plot_path)
    plt.close()

    # q: print to the stdout cumulative fractions for row numbers 100, 500, 1000, 5000, 10000
    print("Cumulative fractions at specified row numbers:")
    for row in [100, 500, 1000, 5000, 10000]:
        if row <= n_rows:
            fraction = report["cumulative_fraction_genus_match"][row - 1]
            print(f"Row {row}: {fraction:.4f}")
        else:
            print(f"Row {row}: out of range")

    

if __name__ == "__main__":
    main()
