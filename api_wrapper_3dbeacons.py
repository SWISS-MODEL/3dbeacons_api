"""Wrapper for the 3D-Beacons API."""

from pathlib import Path
from urllib.parse import urlsplit
import argparse
import json
import re
import sys

import requests


def _parse_args():
    """Get command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)

    metas = {
        "outdir": "<OUTPUT DIRECTORY>",
        "qualifier": "<QUALIFIER>",
    }
    endpoints = ("UniProt", "Ensembl")
    parser.add_argument(
        "endpoint",
        choices=[x.lower() for x in endpoints],
        help="The database to look up entries, available: '"
        + "', '".join(endpoints)
        + "'",
        metavar="<ENDPOINT>",
        type=str.lower,
    )
    # ToDo: [blocker] Turn into file, can we have both? Maybe a JSON-override?
    #       Can I set a flag for that in a workflow?
    parser.add_argument(
        "qualifier", help="ID/ AC to search for", metavar=metas["qualifier"]
    )
    parser.add_argument(
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
    parser.add_argument(
        "--provider",
        choices=providers,
        help="Request from provider only: '"
        + "', '".join(providers)
        + "'. If none is specified, query all providers",
        metavar="<PROVIDER>",
        type=str.lower,
    )
    parser.add_argument(
        "--exclude-provider",
        choices=providers,
        help="Exclude provider: '"
        + "', '".join(providers)
        + "'. If none is specified, query all providers",
        metavar="<PROVIDER>",
        type=str.lower,
    )
    parser.add_argument(
        "--template",
        help="Template PDB ID (4-character code, e.g. '1ake') or SMTL ID with "
        + "assembly and chain name (e.g. 1ake.1.A)",
        metavar="<TEMPLATE>",
        type=str.lower,
    )

    opts = parser.parse_args()

    # check positional arguments
    for req, mta in metas.items():
        value = getattr(opts, req)
        if isinstance(value, str):
            if len(value) == 0:
                parser.error(f"argument of '{mta}' can not be an empty string")
        elif isinstance(value, Path):
            # Empty path points to '.'
            continue
        else:
            parser.error(
                f"Value with unknown type '{type(value).__name__}' found for "
                + f"'{mta}'"
            )
    # check optional arguments
    # check provider & exclude-provider to be different (404es)
    if opts.provider is not None and opts.exclude_provider is not None:
        if opts.provider == opts.exclude_provider:
            parser.error(
                "--provider and --exclude-provider must have different values"
            )
    # validate template format (4-character PDB ID with optional chain/assembly)
    if opts.template is not None:
        # Matches: 1ake, 1ake.1.A, 1ake1A, 1ake_1_A, etc.
        if not re.match(
            r"^[a-z0-9]{4}([._-]?\d+[._-]?[A-Za-z])?$", opts.template
        ):
            parser.error(
                "--template must be a 4-character PDB ID with optional "
                + "chain/assembly (e.g., '1ake', '1ake.1.A', '1ake1A'), got "
                + f"'{opts.template}'"
            )

    return opts


def _assemble_call(qualifier, endpoint, provider, exclude_provider, template):
    """Get the URL for the request assembled."""
    # ToDo: add sequence range
    # ToDo: add CRC64 for uniprot

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


def _main():
    """Run as script."""
    opts = _parse_args()
    call, params = _assemble_call(
        opts.qualifier,
        opts.endpoint,
        opts.provider,
        opts.exclude_provider,
        opts.template,
    )
    response = _send_request(call, params)
    _collect_results(response, opts.endpoint, opts.outdir)

    sys.exit(0)


if __name__ == "__main__":
    _main()

#  LocalWords:  Pylint FASTA
