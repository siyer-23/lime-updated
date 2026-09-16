from __future__ import annotations

import argparse
import csv
import hashlib

from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


# ---------------------------------------------------------------------
# rMATS event definitions
# ---------------------------------------------------------------------

EVENT_TYPES = ("SE", "A3SS", "A5SS", "MXE", "RI")


EVENT_COORD_COLUMNS = {
    "SE": (
        "exonStart_0base",
        "exonEnd",
        "upstreamES",
        "upstreamEE",
        "downstreamES",
        "downstreamEE",
    ),

    "A3SS": (
        "longExonStart_0base",
        "longExonEnd",
        "shortES",
        "shortEE",
        "flankingES",
        "flankingEE",
    ),

    "A5SS": (
        "longExonStart_0base",
        "longExonEnd",
        "shortES",
        "shortEE",
        "flankingES",
        "flankingEE",
    ),

    "MXE": (
        "1stExonStart_0base",
        "1stExonEnd",
        "2ndExonStart_0base",
        "2ndExonEnd",
        "upstreamES",
        "upstreamEE",
        "downstreamES",
        "downstreamEE",
    ),

    "RI": (
        "riExonStart_0base",
        "riExonEnd",
        "upstreamES",
        "upstreamEE",
        "downstreamES",
        "downstreamEE",
    ),
}


INSTANCE_COLUMNS = (
    "universal_event_id",
    "canonical_key",
    "event_type",

    "run_id",
    "run_dir",

    "group1",
    "group2",
    "comparison_label",

    "rmats_event_id",

    "gene_id",
    "gene_symbol",
    "chr",
    "strand",

    "p_value",
    "fdr",
    "delta_psi",
    "is_significant",

    "inc_level_1",
    "inc_level_2",
)


# ---------------------------------------------------------------------
# Dataset metadata
# ---------------------------------------------------------------------

@dataclass(frozen=True)
class RunConfig:
    """
    Metadata describing one rMATS comparison.

    A dataclass is a convenient way to represent a structured object.
    """

    run_id: str
    run_dir: Path

    group1: str
    group2: str

    comparison_label: str


# ---------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------

def dedupe_header(header: list[str]) -> list[str]:
    """
    Make duplicate column names unique.

    Some rMATS output versions contain ID twice in the header.

    Example:

        ID ... ID ...

    becomes:

        ID ... ID__2 ...
    """

    seen: Counter[str] = Counter()

    result: list[str] = []

    for name in header:

        seen[name] += 1

        if seen[name] == 1:
            result.append(name)

        else:
            result.append(f"{name}__{seen[name]}")

    return result


def strip_quotes(value: str) -> str:
    """
    Convert values such as:

        "ENSG000001234"

    into:

        ENSG000001234
    """

    value = value.strip()

    if (
        len(value) >= 2
        and value[0] == '"'
        and value[-1] == '"'
    ):
        return value[1:-1]

    return value


def parse_float(value: str) -> float | None:
    """
    Safely convert an rMATS numeric value into a Python float.

    Missing rMATS values such as NA become None.
    """

    value = value.strip()

    if value in {"", "NA", "NaN", "nan"}:
        return None

    return float(value)


def format_float(value: float | None) -> str:

    if value is None:
        return ""

    return str(value)


# ---------------------------------------------------------------------
# Universal event identity
# ---------------------------------------------------------------------

def canonical_key(
    event_type: str,
    row: dict[str, str],
) -> str:
    """
    Construct the structural identity of an rMATS event.

    Identity is based on:

        event type
        chromosome
        strand
        event-defining coordinates

    It intentionally does NOT use the rMATS ID.
    """

    coord_cols = EVENT_COORD_COLUMNS[event_type]

    required = (
        "chr",
        "strand",
        *coord_cols,
    )

    missing = [
        column
        for column in required
        if column not in row
    ]

    if missing:

        raise ValueError(
            f"{event_type}: missing required columns: "
            f"{', '.join(missing)}"
        )

    parts = [
        event_type,
        strip_quotes(row["chr"]),
        strip_quotes(row["strand"]),
        *(
            strip_quotes(row[column])
            for column in coord_cols
        ),
    ]

    return "|".join(parts)


def universal_event_id(
    event_type: str,
    key: str,
) -> str:
    """
    Convert a potentially long canonical key into a stable short ID.
    """

    digest = hashlib.sha256(
        key.encode("utf-8")
    ).hexdigest()[:16]

    return f"LIME-{event_type}-{digest}"


# ---------------------------------------------------------------------
# Manifest parser
# ---------------------------------------------------------------------

def read_manifest(path: Path) -> list[RunConfig]:

    configs: list[RunConfig] = []

    with path.open(newline="") as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        required = {
            "run_id",
            "run_dir",
            "group1",
            "group2",
            "comparison_label",
        }

        missing = required - set(
            reader.fieldnames or []
        )

        if missing:

            raise ValueError(
                "Manifest missing columns: "
                + ", ".join(sorted(missing))
            )

        for row in reader:

            configs.append(

                RunConfig(
                    run_id=row["run_id"].strip(),

                    run_dir=Path(
                        row["run_dir"].strip()
                    ),

                    group1=row["group1"].strip(),

                    group2=row["group2"].strip(),

                    comparison_label=(
                        row["comparison_label"].strip()
                    ),
                )

            )

    if not configs:
        raise ValueError(
            "Manifest contains no runs."
        )

    run_ids = [
        config.run_id
        for config in configs
    ]

    if len(set(run_ids)) != len(run_ids):

        raise ValueError(
            "run_id values must be unique."
        )

    return configs


# ---------------------------------------------------------------------
# rMATS parser
# ---------------------------------------------------------------------

def iter_rmats_rows(
    path: Path,
) -> Iterable[dict[str, str]]:
    """
    Read an rMATS TSV while validating its structure.

    The row-width check is deliberately strict. A malformed rMATS
    table should cause the program to stop rather than silently
    generating incorrect genomic event IDs.
    """

    with path.open(newline="") as handle:

        reader = csv.reader(
            handle,
            delimiter="\t",
        )

        try:

            raw_header = next(reader)

        except StopIteration:

            raise ValueError(
                f"Empty rMATS file: {path}"
            )

        header = dedupe_header(
            raw_header
        )

        for line_number, fields in enumerate(
            reader,
            start=2,
        ):

            if (
                not fields
                or all(
                    not field.strip()
                    for field in fields
                )
            ):
                continue

            if len(fields) != len(header):

                raise ValueError(
                    f"Malformed tabular row in {path} "
                    f"at line {line_number}: "
                    f"header has {len(header)} columns "
                    f"but row has {len(fields)}. "
                    "Check the source rMATS file."
                )

            yield dict(
                zip(header, fields)
            )


# ---------------------------------------------------------------------
# Significance classification
# ---------------------------------------------------------------------

def significant(
    fdr: float | None,
    delta_psi: float | None,
    fdr_threshold: float,
    min_abs_delta_psi: float,
) -> bool:

    if (
        fdr is None
        or delta_psi is None
    ):
        return False

    return (
        fdr <= fdr_threshold
        and abs(delta_psi)
        >= min_abs_delta_psi
    )


# ---------------------------------------------------------------------
# Load all event instances
# ---------------------------------------------------------------------

def load_instances(
    configs: list[RunConfig],
    count_mode: str,
    fdr_threshold: float,
    min_abs_delta_psi: float,
) -> list[dict[str, str]]:

    instances: list[dict[str, str]] = []

    for config in configs:

        if not config.run_dir.is_dir():

            raise FileNotFoundError(
                "rMATS directory does not exist: "
                f"{config.run_dir}"
            )

        for event_type in EVENT_TYPES:

            path = (
                config.run_dir
                / f"{event_type}.MATS."
                  f"{count_mode}.txt"
            )

            if not path.exists():

                raise FileNotFoundError(
                    "Expected rMATS file not found: "
                    f"{path}"
                )

            for row in iter_rmats_rows(path):

                key = canonical_key(
                    event_type,
                    row,
                )

                uid = universal_event_id(
                    event_type,
                    key,
                )

                p_value = parse_float(
                    row.get(
                        "PValue",
                        "",
                    )
                )

                fdr = parse_float(
                    row.get(
                        "FDR",
                        "",
                    )
                )

                delta_psi = parse_float(
                    row.get(
                        "IncLevelDifference",
                        "",
                    )
                )

                is_sig = significant(
                    fdr=fdr,
                    delta_psi=delta_psi,
                    fdr_threshold=fdr_threshold,
                    min_abs_delta_psi=(
                        min_abs_delta_psi
                    ),
                )

                instances.append({

                    "universal_event_id": uid,

                    "canonical_key": key,

                    "event_type": event_type,

                    "run_id": config.run_id,

                    "run_dir": str(
                        config.run_dir
                    ),

                    "group1": config.group1,

                    "group2": config.group2,

                    "comparison_label": (
                        config.comparison_label
                    ),

                    "rmats_event_id": (
                        strip_quotes(
                            row.get("ID", "")
                        )
                    ),

                    "gene_id": (
                        strip_quotes(
                            row.get("GeneID", "")
                        )
                    ),

                    "gene_symbol": (
                        strip_quotes(
                            row.get(
                                "geneSymbol",
                                "",
                            )
                        )
                    ),

                    "chr": (
                        strip_quotes(
                            row.get("chr", "")
                        )
                    ),

                    "strand": (
                        strip_quotes(
                            row.get(
                                "strand",
                                "",
                            )
                        )
                    ),

                    "p_value": (
                        format_float(p_value)
                    ),

                    "fdr": (
                        format_float(fdr)
                    ),

                    "delta_psi": (
                        format_float(
                            delta_psi
                        )
                    ),

                    "is_significant": (
                        "1"
                        if is_sig
                        else "0"
                    ),

                    "inc_level_1": (
                        row.get(
                            "IncLevel1",
                            "",
                        )
                    ),

                    "inc_level_2": (
                        row.get(
                            "IncLevel2",
                            "",
                        )
                    ),
                })

    return instances


# ---------------------------------------------------------------------
# TSV writer
# ---------------------------------------------------------------------

def write_tsv(
    path: Path,
    rows: list[dict[str, str]],
    columns: Iterable[str],
) -> None:

    columns = list(columns)

    with path.open(
        "w",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=columns,
            extrasaction="ignore",
        )

        writer.writeheader()

        writer.writerows(rows)


# ---------------------------------------------------------------------
# Universal-event catalog
# ---------------------------------------------------------------------

def build_catalog(
    instances: list[dict[str, str]],
    configs: list[RunConfig],
) -> list[dict[str, str]]:

    by_uid: dict[
        str,
        list[dict[str, str]]
    ] = defaultdict(list)

    for row in instances:

        by_uid[
            row["universal_event_id"]
        ].append(row)

    run_ids = [
        config.run_id
        for config in configs
    ]

    catalog: list[
        dict[str, str]
    ] = []

    for uid, rows in sorted(
        by_uid.items()
    ):

        first = rows[0]

        genes = sorted({
            row["gene_id"]
            for row in rows
            if row["gene_id"]
        })

        symbols = sorted({
            row["gene_symbol"]
            for row in rows
            if row["gene_symbol"]
        })

        out = {

            "universal_event_id": uid,

            "canonical_key": (
                first["canonical_key"]
            ),

            "event_type": (
                first["event_type"]
            ),

            "chr": first["chr"],

            "strand": first["strand"],

            "gene_ids": ";".join(
                genes
            ),

            "gene_symbols": ";".join(
                symbols
            ),

            "n_runs_present": str(
                len({
                    row["run_id"]
                    for row in rows
                })
            ),

            "n_runs_significant": str(
                len({
                    row["run_id"]
                    for row in rows
                    if row[
                        "is_significant"
                    ] == "1"
                })
            ),
        }

        for run_id in run_ids:

            run_rows = [
                row
                for row in rows
                if row["run_id"] == run_id
            ]

            out[
                f"present__{run_id}"
            ] = (
                "1"
                if run_rows
                else "0"
            )

            out[
                f"significant__{run_id}"
            ] = (
                "1"
                if any(
                    row[
                        "is_significant"
                    ] == "1"
                    for row in run_rows
                )
                else "0"
            )

            out[
                f"delta_psi__{run_id}"
            ] = (
                run_rows[0]["delta_psi"]
                if run_rows
                else ""
            )

            out[
                f"fdr__{run_id}"
            ] = (
                run_rows[0]["fdr"]
                if run_rows
                else ""
            )

        catalog.append(out)

    return catalog


# ---------------------------------------------------------------------
# Presence / significance pattern annotation
# ---------------------------------------------------------------------

def add_presence_patterns(
    catalog: list[dict[str, str]],
    run_ids: list[str],
) -> None:

    for row in catalog:

        present = [
            run_id
            for run_id in run_ids
            if row[
                f"present__{run_id}"
            ] == "1"
        ]

        significant_runs = [
            run_id
            for run_id in run_ids
            if row[
                f"significant__{run_id}"
            ] == "1"
        ]

        row["presence_pattern"] = (
            "+".join(present)
            if present
            else "none"
        )

        row["significance_pattern"] = (
            "+".join(significant_runs)
            if significant_runs
            else "none"
        )


# ---------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------

def build_overlap_summary(
    catalog: list[dict[str, str]],
) -> list[dict[str, str]]:

    counts = Counter()

    for row in catalog:

        counts[
            (
                row["event_type"],
                row["presence_pattern"],
                row[
                    "significance_pattern"
                ],
            )
        ] += 1

    output = []

    for (
        event_type,
        presence,
        significance_pattern,
    ), n in sorted(
        counts.items()
    ):

        output.append({

            "event_type": event_type,

            "presence_pattern": (
                presence
            ),

            "significance_pattern": (
                significance_pattern
            ),

            "n_events": str(n),
        })

    return output


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def main() -> None:

    parser = argparse.ArgumentParser(

        description=(
            "Build a cross-run catalog "
            "of structurally identical "
            "rMATS events."
        )
    )

    parser.add_argument(
        "--manifest",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--count-mode",
        choices=(
            "JC",
            "JCEC",
        ),
        default="JC",
    )

    parser.add_argument(
        "--fdr",
        type=float,
        default=0.05,
    )

    parser.add_argument(
        "--min-abs-dpsi",
        type=float,
        default=0.0,
    )

    args = parser.parse_args()

    configs = read_manifest(
        args.manifest
    )

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    instances = load_instances(

        configs=configs,

        count_mode=args.count_mode,

        fdr_threshold=args.fdr,

        min_abs_delta_psi=(
            args.min_abs_dpsi
        ),
    )

    run_ids = [
        config.run_id
        for config in configs
    ]

    catalog = build_catalog(
        instances,
        configs,
    )

    add_presence_patterns(
        catalog,
        run_ids,
    )

    overlap = build_overlap_summary(
        catalog
    )

    write_tsv(

        args.output_dir
        / "rmats_event_instances.tsv",

        instances,

        INSTANCE_COLUMNS,
    )

    catalog_columns = [

        "universal_event_id",

        "canonical_key",

        "event_type",

        "chr",

        "strand",

        "gene_ids",

        "gene_symbols",

        "n_runs_present",

        "n_runs_significant",

        "presence_pattern",

        "significance_pattern",
    ]

    for run_id in run_ids:

        catalog_columns.extend([

            f"present__{run_id}",

            f"significant__{run_id}",

            f"delta_psi__{run_id}",

            f"fdr__{run_id}",
        ])

    write_tsv(

        args.output_dir
        / "rmats_event_catalog.tsv",

        catalog,

        catalog_columns,
    )

    write_tsv(

        args.output_dir
        / "rmats_overlap_summary.tsv",

        overlap,

        (
            "event_type",
            "presence_pattern",
            "significance_pattern",
            "n_events",
        ),
    )

    print(
        f"Wrote {len(instances):,} "
        "event instances."
    )

    print(
        f"Wrote {len(catalog):,} "
        "unique structural events."
    )

    print(
        "Output directory: "
        f"{args.output_dir}"
    )


if __name__ == "__main__":
    main()
