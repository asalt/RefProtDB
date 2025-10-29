import csv
from io import StringIO
from pathlib import Path

from collections import Counter, defaultdict

from RefProtDB.cli import cli
from RefProtDB.utils import (
    _fasta_dict_from_file,
    gene_histogram,
    gene_level_fasta,
    group_records_by_gene,
    summarise_gene_groups,
    write_gene_histogram,
    write_gene_summary,
)


def _parse_entries(fasta_text, **kwargs):
    handle = StringIO(fasta_text)
    return list(_fasta_dict_from_file(handle, **kwargs))


FIXTURES = Path(__file__).parent / "data"
FULL_FIXTURE = Path(__file__).parent / "NCBI_refseq_Sgg_tx53354_20251023.fasta"


def _parse_entries_from_path(filename, **kwargs):
    with (FIXTURES / filename).open("r") as handle:
        return list(_fasta_dict_from_file(handle, **kwargs))


def test_specific_header_preserves_full_raw_header_for_wp_accessions():
    fasta = ">WP_12345.1 Hypothetical protein\nMSEQUENCE\n"
    entries = _parse_entries(fasta)
    assert len(entries) == 1
    entry = entries[0]
    assert entry["raw_header"] == ">WP_12345.1 Hypothetical protein"
    assert entry["raw_id"] == "WP_12345.1"
    assert entry["description"] == "Hypothetical protein"
    assert entry["sequence"] == "MSEQUENCE"
    assert entry["protein"] == "WP_12345.1"
    assert entry["geneid"].startswith("meta:")
    assert entry["meta_geneid"] == entry["geneid"]
    assert entry["geneid_source"] == "inferred_name"


def test_specific_header_keeps_full_raw_header_with_uniprot_format():
    fasta = ">sp|Q9BS26|P53_HUMAN Cellular tumor antigen p53\nMPEPTIDE\n"
    entries = _parse_entries(fasta)
    assert len(entries) == 1
    entry = entries[0]
    assert entry["raw_header"] == (
        ">sp|Q9BS26|P53_HUMAN Cellular tumor antigen p53"
    )
    assert entry["raw_id"] == "sp|Q9BS26|P53_HUMAN"
    assert entry["description"] == "Cellular tumor antigen p53"
    assert entry["sp"] == "Q9BS26"
    assert entry["sequence"] == "MPEPTIDE"
    assert entry["protein"] == "sp|Q9BS26|P53_HUMAN"
    assert entry["geneid"].startswith("meta:")
    assert entry["meta_geneid"] == entry["geneid"]
    assert entry["geneid_source"] == "inferred_name"
    assert entry["meta_gene_name"] == "Cellular tumor antigen p53"


def test_generic_header_preserves_full_header_line():
    fasta = ">WP_54321.2 Another description\nMATSEQ\n"
    entries = _parse_entries(fasta, header_search="generic")
    assert len(entries) == 1
    entry = entries[0]
    assert entry["raw_header"] == ">WP_54321.2 Another description"
    assert entry["raw_id"] == "WP_54321.2"
    assert entry["description"] == "Another description"
    assert entry["sequence"] == "MATSEQ"
    assert entry["protein"] == "WP_54321.2"
    assert entry["geneid"].startswith("meta:")


def test_specific_header_real_uniprot_entry_preserves_tokens():
    entries = _parse_entries_from_path("testdata1.fa")
    assert len(entries) == 1
    entry = entries[0]
    expected_header = (FIXTURES / "testdata1.fa").read_text().splitlines()[0]
    assert entry["raw_header"] == expected_header
    assert entry["raw_id"] == expected_header.split()[0][1:]
    assert entry["description"] == ""
    assert entry["ENSP"] == "ENSP00000483084"
    assert entry["ENST"] == "ENST00000611049"
    assert entry["ENSG"] == "ENSG00000152208"
    assert entry["geneid"] == "2895"
    assert entry["taxon"] == "9606"
    assert entry["symbol"] == "GRID2"
    assert entry["sequence"].startswith("ACELMNQGI")
    assert entry["protein"] == entry["raw_id"]
    assert "meta_geneid" not in entry


def test_specific_header_real_wp_entry_keeps_full_description():
    entries = _parse_entries_from_path("testdata2.fa")
    assert len(entries) == 1
    entry = entries[0]
    expected_header = (FIXTURES / "testdata2.fa").read_text().splitlines()[0]
    assert entry["raw_header"] == expected_header
    assert entry["raw_id"] == expected_header.split()[0][1:]
    assert (
        entry["description"]
        == "MULTISPECIES: FtsK/SpoIIIE domain-containing protein [Lactobacillales]"
    )
    assert entry["sequence"].startswith("MLNKGHRIR")
    assert entry["geneid"].startswith("meta:")
    assert entry["protein"] == entry["raw_id"]
    assert entry["geneid_source"] == "inferred_name"


def test_inferred_meta_gene_groups_identical_descriptions():
    entries = _parse_entries_from_path("test_meta_gene.fa")
    assert len(entries) == 6

    def by_description(prefix):
        return {
            e["geneid"]
            for e in entries
            if e["description"].startswith(prefix)
        }

    esa_ids = by_description("type VII secretion protein EsaA")
    effector_ids = by_description("TIGR04197 family type VII secretion effector")
    multi_ids = by_description("MULTISPECIES: TIGR04197 family type VII secretion effector")
    assert len(esa_ids) == 1
    assert len(effector_ids) == 1
    assert len(multi_ids) == 1
    esa_meta = esa_ids.pop()
    effector_meta = effector_ids.pop()
    multi_meta = multi_ids.pop()
    assert esa_meta != effector_meta
    assert effector_meta == multi_meta
    assert esa_meta != multi_meta

    hypothetical = next(
        e for e in entries if e["description"].startswith("hypothetical protein")
    )
    assert hypothetical["geneid"].startswith("meta:p:")
    assert hypothetical["geneid"] == hypothetical["meta_geneid"]
    assert hypothetical["geneid_source"] == "inferred_name"
    assert "meta_gene_name" not in hypothetical
    assert hypothetical["protein"] == hypothetical["raw_id"]

    export_path = FIXTURES / "test_meta_gene.tsv"
    groups = group_records_by_gene(entries)
    summaries = summarise_gene_groups(groups)
    write_gene_summary(summaries, export_path)


def test_meta_gene_inference_on_full_sgg_fixture():
    path = FULL_FIXTURE
    targets = {
        "esaA": "type VII secretion protein EsaA",
        "essB": "type VII secretion protein EssB",
        "tigr": "TIGR04197 family type VII secretion effector",
    }
    groups = defaultdict(list)

    with path.open("r") as handle:
        for entry in _fasta_dict_from_file(handle):
            groups[entry["geneid"]].append(entry)

    def summarise(entries):
        proteins = sorted({e["raw_id"] for e in entries})
        descriptions = sorted({e["description"] for e in entries})
        name = next(
            (e["meta_gene_name"] for e in entries if e.get("meta_gene_name")),
            descriptions[0],
        )
        source = next(
            (e["geneid_source"] for e in entries if e.get("geneid_source")),
            "provided",
        )
        taxons = sorted({e.get("taxon") for e in entries if e.get("taxon")})
        return {
            "geneid": entries[0]["geneid"],
            "meta_gene_name": name,
            "geneid_source": source,
            "count": len(proteins),
            "proteins": "|".join(proteins),
            "descriptions": "|".join(descriptions),
            "taxa": "|".join(taxons),
        }

    summaries = {gid: summarise(entries) for gid, entries in groups.items()}

    def gid_for_prefix(prefix):
        matches = [
            gid
            for gid, entries in groups.items()
            if any(e["description"].startswith(prefix) for e in entries)
        ]
        assert matches, f"No groups found for prefix {prefix}"
        assert len(matches) == 1
        return matches[0]

    esa_gid = gid_for_prefix(targets["esaA"])
    essb_gid = gid_for_prefix(targets["essB"])
    tigr_gid = gid_for_prefix(targets["tigr"])
    assert esa_gid != essb_gid
    assert esa_gid != tigr_gid
    assert essb_gid != tigr_gid

    esa_summary = summaries[esa_gid]
    assert esa_summary["geneid_source"] == "inferred_name"
    assert esa_summary["meta_gene_name"] == targets["esaA"]
    assert esa_summary["count"] >= 10

    essb_summary = summaries[essb_gid]
    assert essb_summary["geneid_source"] == "inferred_name"
    assert essb_summary["meta_gene_name"] == targets["essB"]
    assert essb_summary["count"] >= 5

    tigr_summary = summaries[tigr_gid]
    assert tigr_summary["geneid_source"] == "inferred_name"
    assert tigr_summary["meta_gene_name"] == targets["tigr"]
    assert tigr_summary["count"] >= 10

    multispecies_gid = gid_for_prefix(
        "MULTISPECIES: TIGR04197 family type VII secretion effector"
    )
    assert multispecies_gid == tigr_gid

    export_path = path.with_suffix(".tsv")
    rows = sorted(summaries.values(), key=lambda r: (r["geneid"], r["meta_gene_name"]))
    export_path.parent.mkdir(parents=True, exist_ok=True)
    with export_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            [
                "geneid",
                "meta_gene_name",
                "geneid_source",
                "count",
                "descriptions",
                "proteins",
                "taxa",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    histogram_path = path.with_name(path.name.replace(".fasta", "_hist.tsv"))
    histogram = Counter(summary["count"] for summary in summaries.values())
    with histogram_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            ["proteins_per_gene", "gene_count"],
            delimiter="\t",
        )
        writer.writeheader()
        for size in sorted(histogram):
            writer.writerow(
                {"proteins_per_gene": size, "gene_count": histogram[size]}
            )
