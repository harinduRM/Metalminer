"""CLI wrapper for MetalMiner."""

import argparse
import json

from .core import Config, run_pipeline


def _split_list(values):
    if not values:
        return []
    items = []
    for value in values:
        for part in str(value).split(","):
            part = part.strip()
            if part:
                items.append(part)
    return items


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="MetalMiner CLI")
    parser.add_argument(
        "--target-metals",
        nargs="+",
        default=Config().TARGET_METAL_list,
        help="One or more target metals (space or comma-separated).",
    )
    parser.add_argument(
        "--target-csd-ids",
        nargs="*",
        default=None,
        help="Optional list of CSD IDs to process (space or comma-separated).",
    )
    parser.add_argument("--extraction-method", default=Config().EXTRACTION_METHOD)
    parser.add_argument("--r-factor-limit", type=float, default=Config().R_FACTOR_LIMIT)
    parser.add_argument("--process-limit", type=int, default=Config().PROCESS_LIMIT)
    parser.add_argument(
        "--site-timeout-seconds",
        type=int,
        default=Config().SITE_TIMEOUT_SECONDS,
        help="Optional per-site timeout in seconds. 0 or unset disables timeouts.",
    )
    parser.add_argument(
        "--retry-timeouts",
        action=argparse.BooleanOptionalAction,
        default=Config().RETRY_TIMEOUTS,
        help="Retry timed-out sites when resuming from a spool.",
    )
    parser.add_argument("--visualize", action="store_true", default=Config().VISUALIZE)
    parser.add_argument("--visualization-limit", type=int, default=Config().VISUALIZATION_LIMIT)
    parser.add_argument("--filter-polymeric", action=argparse.BooleanOptionalAction, default=Config().FILTER_POLYMERIC)
    parser.add_argument("--filter-powder", action=argparse.BooleanOptionalAction, default=Config().FILTER_POWDER)
    parser.add_argument("--filter-all-disorder", action=argparse.BooleanOptionalAction, default=Config().FILTER_ALL_DISORDER)
    parser.add_argument("--filter-primary-disorder", action=argparse.BooleanOptionalAction, default=Config().FILTER_PRIMARY_DISORDER)
    parser.add_argument("--correct-primary-disorder", action=argparse.BooleanOptionalAction, default=Config().CORRECT_PRIMARY_DISORDER)
    parser.add_argument("--hydrogen-addition-method", default=Config().Hydrogen_Addition_method)
    parser.add_argument("--disorder-resolve-method", default=Config().DISORDER_RESOLVE_METHOD)
    parser.add_argument("--num-metal-layers", type=int, default=Config().num_metal_layers)
    parser.add_argument(
        "--geometric-radius",
        type=float,
        default=Config().Geometric_radius,
        help="Radius used for geometric extraction (only when extraction method is Geometric).",
    )
    parser.add_argument(
        "--oxidation-states-filter",
        nargs="+",
        default=Config().OXIDATION_STATES_FILTER,
        help="List of allowed oxidation states or 'all'.",
    )
    parser.add_argument("--get-abstract", action=argparse.BooleanOptionalAction, default=Config().GET_ABSTRACT)
    parser.add_argument("--edit-manual", action=argparse.BooleanOptionalAction, default=Config().Edit_manual)
    parser.add_argument("--limit-nm-nm", type=float, default=Config().Extraction_Cut_off_distances["LIMIT_NM_NM"])
    parser.add_argument("--limit-m-nm", type=float, default=Config().Extraction_Cut_off_distances["LIMIT_M_NM"])
    parser.add_argument("--limit-h-x", type=float, default=Config().Extraction_Cut_off_distances["LIMIT_H_X"])
    parser.add_argument(
        "--extraction-cutoff-json",
        default=None,
        help="Optional JSON string to override extraction cut-off distances.",
    )
    parser.add_argument(
        "--metalloligands",
        nargs="+",
        default=Config().metalloligands,
        help="List of metalloligand elements (space or comma-separated).",
    )
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    target_metals = _split_list(args.target_metals)
    target_csd_ids = _split_list(args.target_csd_ids) if args.target_csd_ids else None
    oxidation_states = _split_list(args.oxidation_states_filter)
    metalloligands = _split_list(args.metalloligands)

    extraction_cutoffs = {
        "LIMIT_NM_NM": args.limit_nm_nm,
        "LIMIT_M_NM": args.limit_m_nm,
        "LIMIT_H_X": args.limit_h_x,
    }
    if args.extraction_cutoff_json:
        extraction_cutoffs = json.loads(args.extraction_cutoff_json)

    config = Config(
        TARGET_METAL_list=target_metals,
        target_csd_ids=target_csd_ids,
        EXTRACTION_METHOD=args.extraction_method,
        R_FACTOR_LIMIT=args.r_factor_limit,
        PROCESS_LIMIT=args.process_limit,
        SITE_TIMEOUT_SECONDS=args.site_timeout_seconds,
        RETRY_TIMEOUTS=args.retry_timeouts,
        VISUALIZE=args.visualize,
        VISUALIZATION_LIMIT=args.visualization_limit,
        FILTER_POLYMERIC=args.filter_polymeric,
        FILTER_POWDER=args.filter_powder,
        FILTER_ALL_DISORDER=args.filter_all_disorder,
        FILTER_PRIMARY_DISORDER=args.filter_primary_disorder,
        CORRECT_PRIMARY_DISORDER=args.correct_primary_disorder,
        Hydrogen_Addition_method=args.hydrogen_addition_method,
        DISORDER_RESOLVE_METHOD=args.disorder_resolve_method,
        num_metal_layers=args.num_metal_layers,
        OXIDATION_STATES_FILTER=oxidation_states,
        GET_ABSTRACT=args.get_abstract,
        Edit_manual=args.edit_manual,
        Geometric_radius=args.geometric_radius,
        Extraction_Cut_off_distances=extraction_cutoffs,
        metalloligands=metalloligands,
    )
    return run_pipeline(config)


if __name__ == "__main__":
    main()
