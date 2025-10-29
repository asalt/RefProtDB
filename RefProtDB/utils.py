from __future__ import print_function

import os
import re
from datetime import datetime
from functools import wraps
import difflib
import textwrap
import logging
import hashlib
import csv
from collections import defaultdict
from pathlib import Path

import six

import numpy as np
import pandas as pd
import click

from .containers import FASTA_FMT

logger = logging.getLogger(__name__)

end = "..."

header_pat = re.compile(r"(\w+)\|([\w-]+)\|([\w\-]+)?")


_GENERIC_DESCRIPTIONS = {
    "hypothetical protein",
    "putative protein",
    "uncharacterized protein",
    "predicted protein",
}


def _strip_trailing_organism(description):
    if not description:
        return ""
    return re.sub(r"\s*\[[^\]]+\]\s*$", "", description).strip()


def _clean_description(description):
    desc = _strip_trailing_organism(description)
    if not desc:
        return ""
    if desc.upper().startswith("MULTISPECIES:"):
        desc = desc.split(":", 1)[1].strip()
    desc = re.sub(r"(,?\s+partial)$", "", desc, flags=re.IGNORECASE).strip()
    return re.sub(r"\s+", " ", desc).strip()


def _slugify(text, max_length=64):
    slug = re.sub(r"[^a-z0-9]+", "_", text.lower())
    slug = re.sub(r"_+", "_", slug).strip("_")
    if not slug:
        return ""
    return slug[:max_length]


def _meta_gene_key(description, raw_id):
    raw_id = raw_id or ""
    cleaned = _clean_description(description)
    if not cleaned:
        return f"raw:{raw_id}", raw_id
    normalized = cleaned.lower()
    if normalized in _GENERIC_DESCRIPTIONS:
        return f"raw:{raw_id}", raw_id
    return f"name:{normalized}", cleaned


def _meta_gene_identifier(key, label):
    label = label or key
    prefix = "n" if key.startswith("name:") else "p"
    slug = _slugify(label)
    if not slug:
        slug = hashlib.md5(label.encode("utf-8")).hexdigest()[:12]
    return f"meta:{prefix}:{slug}"


def _best_description(entries):
    for entry in entries:
        if entry.get("meta_gene_name"):
            return entry["meta_gene_name"]
    for entry in entries:
        if entry.get("description"):
            return entry["description"]
    return ""


def group_records_by_gene(records):
    groups = defaultdict(list)
    for record in records:
        groups[record["geneid"]].append(record)
    return groups


def summarise_gene_group(geneid, entries):
    description = _best_description(entries)
    proteins = sorted(
        {
            entry.get("raw_id") or entry.get("ref")
            for entry in entries
            if entry.get("raw_id") or entry.get("ref")
        }
    )
    descriptions = sorted({entry.get("description") for entry in entries if entry.get("description")})
    source = next((entry.get("geneid_source") for entry in entries if entry.get("geneid_source")), "provided")
    taxa = sorted({entry.get("taxon") for entry in entries if entry.get("taxon")})
    return {
        "geneid": geneid,
        "meta_gene_name": description,
        "geneid_source": source,
        "count": len(proteins),
        "descriptions": "|".join(descriptions),
        "proteins": "|".join(proteins),
        "taxa": "|".join(taxa),
    }


def summarise_gene_groups(groups):
    return {geneid: summarise_gene_group(geneid, entries) for geneid, entries in groups.items()}


def write_gene_summary(summaries, path):
    if not summaries:
        return
    rows = sorted(summaries.values(), key=lambda r: (r["geneid"], r["meta_gene_name"]))
    fieldnames = [
        "geneid",
        "meta_gene_name",
        "geneid_source",
        "count",
        "descriptions",
        "proteins",
        "taxa",
    ]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def gene_histogram(summaries):
    histogram = defaultdict(int)
    for summary in summaries.values():
        histogram[summary["count"]] += 1
    return dict(sorted(histogram.items()))


def write_gene_histogram(histogram, path):
    if not histogram:
        return
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["proteins_per_gene", "gene_count"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for size, count in histogram.items():
            writer.writerow({"proteins_per_gene": size, "gene_count": count})


def _aggregate_gene_metadata(geneid, entries):
    description = _best_description(entries)
    proteins = sorted(
        {
            entry.get("raw_id") or entry.get("ref")
            for entry in entries
            if entry.get("raw_id") or entry.get("ref")
        }
    )
    refs = sorted(
        {
            entry.get("ref") or entry.get("raw_id")
            for entry in entries
            if entry.get("ref") or entry.get("raw_id")
        }
    )
    gis = sorted({entry.get("gi") for entry in entries if entry.get("gi")})
    symbols = sorted({entry.get("symbol") for entry in entries if entry.get("symbol")})
    taxa = sorted({entry.get("taxon") for entry in entries if entry.get("taxon")})
    homologenes = sorted({entry.get("homologene") for entry in entries if entry.get("homologene")})
    representative = max(entries, key=lambda entry: len(entry.get("sequence", "")))
    sequence = representative.get("sequence", "")
    return {
        "geneid": geneid,
        "description": description,
        "sequence": sequence,
        "ref": "|".join(refs) if refs else representative.get("raw_id", ""),
        "gi": "|".join(gis),
        "symbol": "|".join(symbols),
        "taxon": "|".join(taxa),
        "homologene": "|".join(homologenes),
        "proteins": "|".join(proteins),
    }


def gene_level_fasta(
    fasta_path,
    output_path,
    header_search="specific",
    summary_path=None,
    histogram_path=None,
):
    records = list(fasta_dict_from_file(fasta_path, header_search=header_search))
    groups = group_records_by_gene(records)
    summaries = summarise_gene_groups(groups)
    aggregated = {
        geneid: _aggregate_gene_metadata(geneid, entries)
        for geneid, entries in groups.items()
    }
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as handle:
        for geneid in sorted(aggregated):
            meta = aggregated[geneid]
            seq = "\n".join(textwrap.wrap(meta["sequence"], width=70))
            record = FASTA_FMT.format(
                ref=meta["ref"],
                gis=meta["gi"],
                gid=meta["geneid"],
                homologene=meta["homologene"],
                taxons=meta["taxon"],
                description=meta["description"],
                seq=seq,
                symbols=meta["symbol"],
            )
            handle.write(record)
            handle.write("\n")
    if summary_path:
        write_gene_summary(summaries, summary_path)
    if histogram_path:
        histogram = gene_histogram(summaries)
        write_gene_histogram(histogram, histogram_path)
    return aggregated, summaries


def print_msg(*msg):
    def deco(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            if six.PY3:
                print(datetime.now(), ":", *msg, end=end, flush=True)
            elif six.PY2:
                print(datetime.now(), ":", *msg, end=end)
            result = func(*args, **kwargs)
            print("done.")
            return result

        return wrapper

    return deco


def parse_header(header):
    """Parses a FASTA header, extracting relevant identifiers."""
    TOKEEP = {
        "gi",
        "ref",
        "geneid",
        "symbol",
        "taxon",
        "sp",
        "spid",
        "spacc",
        "spref",
        "ENSG",
        "ENST",
        "ENSP",
    }
    keys = header_pat.findall(header)
    header_data = {"raw_header": header}
    uniprot_id = None

    for key, value, extra in keys:
        if key.startswith("rev"):
            continue  # Ignore reverse/decoy entries
        if key in {"sp", "tr"}:  # UniProt identifiers
            uniprot_id = value
            header_data["uniprot_id"] = uniprot_id
            if extra:
                header_data["uniprot_name"] = extra
        elif key in TOKEEP:
            header_data[key] = value

    # Fetch Entrez ID if only UniProt ID is present
    if uniprot_id and "geneid" not in header_data:
        entrez_id = fetch_entrez_from_uniprot(uniprot_id)
        if entrez_id:
            header_data["geneid"] = entrez_id

    return header_data


# heavily inspired and borrowed
# by https://github.com/jorvis/biocode/blob/master/lib/biocode/utils.py#L149
# Joshua Orvis


def _fasta_dict_from_file(
    file_object,
    header_search="specific",
    decoy_prefix="rev",
    remove_decoys=True,
):
    """
    Reads a file of FASTA entries and returns a dict for each entry.
    It parsers the headers such that `>gi|12345|ref|NP_XXXXX| Description`
    returns as `{gi: 12345, ref: 'NP_XXXXX', description : 'Description'}`
    The sequence has all whitespace removed.

    :header_search: One of `specific` or `generic`
    if `specific` tries to parse the header
    if `general` just returns each whitespace separated header
    """

    current_id = dict()
    current_seq = ""
    current_header = None
    pat = re.compile(r">(\S+)\s*(.*)")
    # old:
    #header_pat = re.compile(r"(\w+)\|([\w\.-]+)?")
    # new allows : inside a header name 
    header_pat = re.compile(r"(\w+)\|([\w\.\-\:]+)?")
    # Explanation:
    # (\w+)       : Capture one or more word characters (letters, digits, or underscores).
    # \|          : Match the pipe character '|'.
    # ([\w\.-\:]+)? : Capture zero or more word characters, dots, or dashes or colons. This part is optional.

    inferred_geneids = dict()

    def assign_meta_geneid(record):
        if record.get("geneid"):
            return
        raw_id_value = record.get("raw_id") or ""
        key, label = _meta_gene_key(record.get("description", ""), raw_id_value)
        meta_geneid = inferred_geneids.get(key)
        if not meta_geneid:
            meta_geneid = _meta_gene_identifier(key, label)
            inferred_geneids[key] = meta_geneid
        record["geneid"] = meta_geneid
        record["meta_geneid"] = meta_geneid
        record.setdefault("geneid_source", "inferred_name")
        if key.startswith("name:"):
            record.setdefault("meta_gene_name", label)

    def parse_header(header, raw_header, raw_id, pairs=True):
        TOKEEP = [
            "gi",
            "ref",
            "geneid",
            "symbol",
            "taxon",
            "sp",
            "spid",
            "spacc",
            "spref",
            "ENSG",
            "ENST",
            "ENSP",
        ]  # no keys other than these are kept. if starts with rev, also dropped
        # if "P00330" in header:
        #     import ipdb

        #     ipdb.set_trace()
        keys = header_pat.findall(header)
        header_data = dict(raw_header=raw_header, raw_id=raw_id)
        invalid_keys = set()
        for key, value in keys:
            # if key[1] == "geneid":
            if key.startswith(decoy_prefix) and remove_decoys:
                # decoy
                continue
            if key not in TOKEEP:
                if key not in invalid_keys:  # only warn once
                    invalid_keys.add(key)
                    # logger.warning(f"key {key} not in {TOKEEP}")
                    continue
                # You can continue to the next iteration if the key is not in TOKEEP
                # skipping
            elif key in TOKEEP:
                header_data[key] = value

        # gi -> ProteinGI #, ref -> NP_XXXX
        return header_data

    for line in file_object:
        line = line.rstrip()
        m = pat.search(line)
        if m:
            ## new residue line matched, purge the existing one, if not the first
            if current_id:
                ## remove all whitespace and save
                current_seq = "".join(current_seq.split())
                current_id["sequence"] = current_seq
                yield current_id
                # current_id.clear()  # this is actually bad for list comprehensions
                # as it returns empty dictionaries

            current_seq = ""
            header = m.group(1)
            raw_header = line
            raw_id = m.group(1)
            if header_search == "specific":
                current_id = parse_header(header, raw_header, raw_id)
            elif header_search == "generic":
                current_id = dict(raw_header=raw_header, raw_id=raw_id)
            current_id["description"] = m.group(2)
            current_id.setdefault("protein", current_id.get("ref") or raw_id)
            assign_meta_geneid(current_id)

        else:
            ## python 2.6+ makes string concatenation amortized O(n)
            ##  http://stackoverflow.com/a/4435752/1368079
            current_seq += str(line)

    ## don't forget the last one
    current_seq = "".join(current_seq.split())
    current_id["sequence"] = current_seq
    yield current_id


def fasta_dict_from_file(fasta_file, header_search="specific"):
    with open(fasta_file, "r") as f:
        # yield from _fasta_dict_from_file(f, header_search=header_search)
        for v in _fasta_dict_from_file(f, header_search=header_search):
            yield v


fasta_dict_from_file.__doc__ = _fasta_dict_from_file.__doc__


def convert_tab_to_fasta(
    tabfile,
    geneid=None,
    ref=None,
    gi=None,
    homologene=None,
    taxon=None,
    description=None,
    sequence=None,
    symbol=None,
):
    HEADERS = {
        "geneid": geneid,
        "ref": ref,
        "gi": gi,
        "homologene": homologene,
        "taxon": taxon,
        "description": description,
        "sequence": sequence,
        "symbol": symbol,
    }
    CUTOFF = 0.35
    df = pd.read_table(tabfile, dtype=str)
    choices = click.Choice([x for y in [df.columns, ["SKIP"]] for x in y])
    col_names = dict()
    for h, v in HEADERS.items():
        if v is not None or v == "SKIP":
            col_names[h] = v
            continue

        closest_match = difflib.get_close_matches(h, df.columns, n=1, cutoff=CUTOFF)
        if closest_match and h != "homologene":
            col_names[h] = closest_match[0]
        else:
            print("Can not find header match for :", h)
            print("Choose from :", " ".join(choices.choices))
            choice = click.prompt(
                "Selected the correct name or SKIP",
                type=choices,
                show_default=True,
                err=True,
                default="SKIP",
            )
            col_names[h] = choice
            print()
    for _, series in df.iterrows():
        row = series.to_dict()
        gid = row.get(col_names["geneid"], "")
        ref = row.get(col_names["ref"], "")
        hid = row.get(col_names["homologene"], "")
        gi = row.get(col_names["gi"], "")
        taxon = row.get(col_names["taxon"], "")
        desc = row.get(col_names["description"], " ")
        seq = "\n".join(textwrap.wrap(row.get(col_names["sequence"], ""), width=70))
        symbol = row.get(col_names["symbol"], "")

        hid = int(hid) if hid and not np.isnan(hid) else ""
        try:
            gid = int(gid)
        except ValueError:
            pass

        r = dict(
            gid=gid,
            ref=ref,
            taxons=taxon,
            gis=gi,
            homologene=hid,
            description=desc,
            seq=seq,
            symbols=symbol,
        )
        yield (FASTA_FMT.format(**r))


def sniff_fasta(fasta):
    nrecs = 1
    fasta = fasta_dict_from_file(fasta)
    counter = 0
    REQUIRED = ("ref", "sequence")
    while counter < nrecs:
        for rec in fasta:
            if any(x not in rec for x in REQUIRED):
                raise ValueError("Invalid input FASTA")
        counter += 1
    return 0  # all good
