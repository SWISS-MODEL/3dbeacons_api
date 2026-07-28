"""Wrapper for the 3D-Beacons API."""

from pathlib import Path
from urllib.parse import urlsplit
import argparse
import csv
import json
import re
import sys

import requests


def _parse_args():
    """Get command line arguments."""
    # pylint: disable=too-many-branches
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="mode", required=True)

    # Single query mode
    single_parser = subparsers.add_parser("single", help="Query a single entry")
    metas = {
        "outdir": "<OUTPUT DIRECTORY>",
        "qualifier": "<QUALIFIER>",
        "endpoint": "<ENDPOINT>",
    }
    endpoints = ("UniProt", "Ensembl")
    single_parser.add_argument(
        "endpoint",
        choices=[x.lower() for x in endpoints],
        help="The database to look up entries, available: '"
        + "', '".join(endpoints)
        + "'",
        metavar="<ENDPOINT>",
        type=str.lower,
    )
    single_parser.add_argument(
        "qualifier",
        help="ID/AC to search for",
        metavar=metas["qualifier"],
    )
    single_parser.add_argument(
        "outdir",
        help="Directory to store results in",
        metavar=metas["outdir"],
        type=Path,
    )
    providers = [
        "alphafill",
        "alphafold",
        "hegelab",
        "isoformio",
        "levylab",
        "modelarchive",
        "pdbe",
        "ped",
        "sasbdb",
        "swissmodel",
    ]
    single_parser.add_argument(
        "--provider",
        choices=providers,
        help="Request from provider only: '"
        + "', '".join(providers)
        + "'. If none is specified, query all providers",
        metavar="<PROVIDER>",
        type=str.lower,
    )
    single_parser.add_argument(
        "--exclude-provider",
        choices=providers,
        help="Exclude provider: '"
        + "', '".join(providers)
        + "'. If none is specified, query all providers",
        metavar="<PROVIDER>",
        type=str.lower,
    )
    single_parser.add_argument(
        "--template",
        help="Template PDB ID (4-character code, e.g. '1ake') or SMTL ID with "
        + "assembly and chain name (e.g. 1ake.1.A)",
        metavar="<TEMPLATE>",
        type=str.lower,
    )
    single_parser.add_argument(
        "--sequence-range",
        help="Specify a UniProt sequence residue range (format: 'start-end', "
        + "e.g., '1-100')",
        metavar="<RANGE>",
    )
    single_parser.add_argument(
        "--uniprot-checksum",
        help="CRC64 checksum of the UniProt sequence (16-character hex value)",
        metavar="<CHECKSUM>",
        type=str.upper,
    )

    # Batch mode
    batch_parser = subparsers.add_parser(
        "batch", help="Process multiple queries from a TSV file"
    )
    batch_parser.add_argument(
        "tsv_file",
        help="Path to TSV query file with columns: endpoint, qualifier, "
        + "and optional parameters (provider, exclude_provider, template, "
        + "sequence_range, uniprot_checksum)",
        metavar="<TSV_FILE>",
        type=Path,
    )
    batch_parser.add_argument(
        "outdir",
        help="Directory to store results in",
        metavar="<OUTPUT_DIRECTORY>",
        type=Path,
    )

    opts = parser.parse_args()

    # Validate single-query mode
    if opts.mode == "single":
        # check provider & exclude-provider to be different (404es)
        if opts.provider is not None and opts.exclude_provider is not None:
            if opts.provider == opts.exclude_provider:
                single_parser.error(
                    "--provider and --exclude-provider must have different "
                    + "values"
                )
        # validate template format (PDB ID with optional chain/assembly)
        if opts.template is not None:
            # Matches: 1ake, 1ake.1.A, 1ake1A, 1ake_1_A, etc.
            if not re.match(
                r"^[a-z0-9]{4}([._-]?\d+[._-]?[A-Za-z])?$", opts.template
            ):
                single_parser.error(
                    "--template must be a 4-character PDB ID with optional "
                    + "chain/assembly (e.g., '1ake', '1ake.1.A', '1ake1A'), "
                    + f"got '{opts.template}'"
                )
        # validate sequence range format (start-end)
        if opts.sequence_range is not None:
            if not re.match(r"^\d+-\d+$", opts.sequence_range):
                single_parser.error(
                    "--sequence-range must be in format 'start-end' (e.g., "
                    + f"'1-100'), got '{opts.sequence_range}'"
                )
        # validate CRC64 checksum format (16 hex characters)
        if opts.uniprot_checksum is not None:
            if not re.match(r"^[A-F0-9]{16}$", opts.uniprot_checksum):
                single_parser.error(
                    "--uniprot-checksum must be a 16-character CRC64 hex value "
                    + "(e.g., '1A2B3C4D5E6F7890'), got "
                    + f"'{opts.uniprot_checksum}'"
                )
    # Validate batch mode
    elif opts.mode == "batch":
        if not opts.tsv_file.is_file():
            batch_parser.error(f"Query file not found: '{opts.tsv_file}'")

    return opts


def _parse_query_file(filepath):
    """Parse TSV query file with flexible column order (streaming).

    Args:
        filepath: Path to TSV query file

    Yields:
        Tuples of (endpoint, qualifier, params_dict)

    Raises:
        ValueError: If file format is invalid
    """

    def _filtered_line(file_handle):
        """Generator that yields non-comment, non-empty lines."""
        for line in file_handle:
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                yield stripped

    try:
        with open(filepath, encoding="ascii") as fh:
            reader = csv.DictReader(_filtered_line(fh), delimiter="\t")

            # Validate header before processing rows
            required = {"endpoint", "qualifier"}
            missing = required - set(reader.fieldnames)
            if missing:
                raise ValueError(
                    f"Query file must have {required} columns in header. "
                    f"Found: {', '.join(reader.fieldnames)}"
                )

            # Process data rows
            for row_num, row in enumerate(reader, start=1):
                endpoint = row.get("endpoint", "").strip().lower()
                if not endpoint:
                    raise ValueError(f"Row {row_num}: empty endpoint value")
                if endpoint not in ("uniprot", "ensembl"):
                    raise ValueError(
                        f"Row {row_num}: invalid endpoint '{endpoint}'. "
                        + "Must be 'uniprot' or 'ensembl'"
                    )

                qualifier = row.get("qualifier", "").strip()
                if not qualifier:
                    raise ValueError(f"Row {row_num}: empty qualifier value")

                params = {
                    key: value.strip()
                    for key, value in row.items()
                    if key not in ("endpoint", "qualifier")
                    and value
                    and value.strip()
                }

                yield (endpoint, qualifier, params)

    except FileNotFoundError as err:
        raise ValueError(f"Query file not found: '{filepath}'") from err
    except csv.Error as err:
        raise ValueError(f"Invalid TSV format in '{filepath}': {err}") from err


def _assemble_call(
    qualifier,
    endpoint,
    provider,
    exclude_provider,
    template,
    sequence_range,
    uniprot_checksum,
):
    """Get the URL for the request assembled."""
    # pylint: disable=too-many-positional-arguments,too-many-arguments
    url = (
        f"https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/v2/{endpoint}"
        + f"/summary/{qualifier}.json"
    )
    params = {}
    if provider is not None:
        params["provider"] = provider
    if exclude_provider is not None:
        params["exclude_provider"] = exclude_provider
    if template is not None:
        params["template"] = template
    if sequence_range is not None:
        params["range"] = sequence_range
    if uniprot_checksum is not None:
        params["uniprot_checksum"] = uniprot_checksum

    return url, params


def _send_request(call, params):
    """Send a request to the 3D-Beacons API and wait for the response."""
    response = requests.get(
        call,
        headers={"Accept": "application/json"},
        params=params,
        timeout=360,
    )
    response.raise_for_status()

    return response.json()


_TRANSLATE_TABLE = str.maketrans(".:", "__")


def _make_coordfile_name(mdl_info):
    """Assemble a name for a coordinate file."""
    # some providers deliver non-unique file name, use model identifier instead
    if mdl_info["provider"].upper() in ("SWISS-MODEL",):
        mdl_id = mdl_info["model_identifier"].translate(_TRANSLATE_TABLE)
    else:
        mdl_id = str(Path(urlsplit(mdl_info["model_url"]).path).name)
    fname = Path(
        f"{mdl_info['provider']}__{mdl_info['model_format'].upper()}"
        + f"__{mdl_id}"
    )
    # make sure file has an extension
    if not fname.suffix:
        if mdl_info["model_format"].upper() == "MMCIF":
            fname = fname.with_suffix(".cif")
        else:
            raise NotImplementedError(
                f"Unknown model format '{mdl_info['model_format']}' found, "
                + "needs to be added to the wrapper script."
            )

    return fname, mdl_id


def _fetch_file(mdl_info, outdir):
    response = requests.get(mdl_info["model_url"], timeout=360)
    response.raise_for_status()
    fname, mid = _make_coordfile_name(mdl_info)
    with open(outdir / fname, "wb") as mfh:
        for chunk in response.iter_content(chunk_size=8192):
            mfh.write(chunk)
    return Path(mid).stem


def _collect_results(response, endpoint, output_dir):
    """Get JSON data and coordinate files."""
    # write entry info as JSON
    out_path = output_dir / "INFO"
    out_path.mkdir(parents=True, exist_ok=True)
    with open(
        out_path / f"{endpoint}_{response['uniprot_entry']['ac']}.json",
        "w",
        encoding="utf8",
    ) as jfh:
        json.dump(response["uniprot_entry"], jfh)
    # get coordinate files & per file a summary JSON
    for mdl in response["structures"]:
        mdl = mdl["summary"]
        # coordinate file
        out_path = output_dir / "MODELS"
        out_path.mkdir(parents=True, exist_ok=True)
        mid = _fetch_file(mdl, out_path)
        # summary file in SUMMARIES/provider/
        summaries_path = output_dir / "SUMMARIES"
        summaries_path.mkdir(parents=True, exist_ok=True)
        with open(
            summaries_path / f"{mdl['provider']}__{mid}.json",
            "w",
            encoding="utf8",
        ) as jfh:
            json.dump(mdl, jfh)


def _run_query(
    qualifier,
    endpoint,
    out_dir,
    provider,
    exclude_provider,
    template,
    sequence_range,
    uniprot_checksum,
):
    """Try to fetch entry data"""
    # pylint: disable=too-many-positional-arguments,too-many-arguments
    call, api_params = _assemble_call(
        qualifier,
        endpoint,
        provider,
        exclude_provider,
        template,
        sequence_range,
        uniprot_checksum,
    )
    try:
        response = _send_request(call, api_params)
    except requests.exceptions.HTTPError:  # Entry not found
        pass
    else:
        _collect_results(response, endpoint, out_dir)


def _main():
    """Run as script."""
    opts = _parse_args()

    if opts.mode == "batch":
        # Batch mode: process TSV file (streaming, memory efficient)
        # Process each query
        queries_processed = False
        for endpoint, qualifier, params in _parse_query_file(opts.tsv_file):
            queries_processed = True
            qual_outdir = opts.outdir / qualifier
            _run_query(
                qualifier,
                endpoint,
                qual_outdir,
                params.get("provider"),
                params.get("exclude_provider"),
                params.get("template"),
                params.get("sequence_range"),
                params.get("uniprot_checksum"),
            )

        if not queries_processed:
            print("Error: Query file contains no data rows", file=sys.stderr)
            sys.exit(1)
    else:
        # Single query mode
        _run_query(
            opts.qualifier,
            opts.endpoint,
            opts.outdir,
            opts.provider,
            opts.exclude_provider,
            opts.template,
            opts.sequence_range,
            opts.uniprot_checksum,
        )

    sys.exit(0)


if __name__ == "__main__":
    _main()

#  LocalWords:  Pylint FASTA
