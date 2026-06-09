"""Wrapper for the 3D-Beacons API."""

from pathlib import Path
from urllib.parse import urlsplit
import argparse
import json
import sys

import requests


def _parse_args():
    """Get command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)

    metas = {
        "outdir": "<OUTPUT DIRECTORY>",
        "qualifier": "<QUALIFIER>",
    }
    endpoints = ("UniProt", "Ensemble")
    parser.add_argument(
        "endpoint",
        choices=[x.lower() for x in endpoints],
        help="The database to look up entries, available: '"
        + "','".join(endpoints)
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
            raise RuntimeError(
                f"Value with unknown type '{type(value).__name__}' found for "
                + f"'{mta}'"
            )

    return opts


def _assemble_call(qualifier, endpoint):
    """Get the URL for the request assembled."""
    # ToDo: Extend with provider
    # ToDo: add template
    # ToDo: add sequence range
    # ToDo: add provider to exclude
    # ToDo: add CRC64 for uniprot
    url = (
        f"https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/v2/{endpoint}"
        + f"/summary/{qualifier}.json"
    )

    return url


def _send_request(call):
    """Send a request to the 3D-Beacons API and wait for the response."""
    response = requests.get(
        call,
        headers={"Accept": "application/json"},
        timeout=360,
    )
    response.raise_for_status()

    return response.json()


_TRANSLATE_TABLE = str.maketrans(".:", "__")


def _make_coordfile_name(mdl_info):
    """Assemble a name for a coordinate file."""
    # some providers deliver non-unique file name, use model identifier instead
    if mdl_info["provider"].upper() in ("SWISS-MODEL",):
        fname = Path(mdl_info["model_identifier"].translate(_TRANSLATE_TABLE))
    else:
        fname = Path(Path(urlsplit(mdl_info["model_url"]).path).name)
    # make sure file has an extension
    if not fname.suffix:
        if mdl_info["model_format"].upper() == "MMCIF":
            fname = fname.with_suffix(".cif")
        else:
            raise NotImplementedError(
                f"Unknown model format '{mdl_info['model_format']}' found, "
                + "needs to be added to the wrapper script."
            )

    return fname


def _fetch_file(mdl_info, outdir):
    response = requests.get(mdl_info["model_url"], timeout=360)
    response.raise_for_status()
    fname = _make_coordfile_name(mdl_info)
    with open(outdir / fname, "wb") as mfh:
        for chunk in response.iter_content(chunk_size=8192):
            mfh.write(chunk)
    return fname.stem


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
        # in provider dir
        out_path = output_dir / mdl["provider"]
        out_path.mkdir(parents=True, exist_ok=True)
        # coordinate file
        fname = _fetch_file(mdl, out_path)
        # summary file
        with open(
            out_path / f"{fname}_summary.json",
            "w",
            encoding="utf8",
        ) as jfh:
            json.dump(mdl, jfh)


def _main():
    """Run as script."""
    opts = _parse_args()
    call = _assemble_call(opts.qualifier, opts.endpoint)
    response = _send_request(call)
    _collect_results(response, opts.endpoint, opts.outdir)

    sys.exit(0)


if __name__ == "__main__":
    _main()

#  LocalWords:  Pylint FASTA
