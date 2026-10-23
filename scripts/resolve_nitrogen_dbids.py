#!/usr/bin/env python3
"""Resolve supplementary table 6 DBids against NCBI and KEGG.

This script is deliberately conservative: a taxon prefix alone does not resolve
the associated locus tag, and opaque identifiers remain unresolved unless an
authoritative namespace mapping is available.
"""

from __future__ import annotations

import csv
import json
import re
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

from openpyxl import load_workbook


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "references/lv_2026/table 6.xlsx"
OUTPUT = ROOT / "data/processed/nitrogen_dbid_mapping.tsv"
AUDIT = ROOT / "docs/nitrogen_mapping_audit.md"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
USER_AGENT = "PRJNA1399466-reproducibility-audit/1.0"
HEX_RE = re.compile(r"^[0-9a-f]{32}$")
INTEGER_RE = re.compile(r"^\d+$")
COMPOSITE_RE = re.compile(r"^(\d+)\.(.+)$")

FIGURE7_CATEGORIES = [
    "Anammox",
    "Nitrification",
    "Auxiliary",
    "Organic degradation and synthesis",
    "Dissimilatory nitrate reduction",
    "Assimilatory nitrate reduction",
    "Denitrification",
    "Denitrification/Dissimilatory nitrate reduction",
    "Nitrogen fixation",
]


def fetch(url: str, *, attempts: int = 4) -> bytes:
    """Fetch a public record politely, retrying transient failures."""
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(request, timeout=45) as response:
                result = response.read()
            time.sleep(0.34)
            return result
        except Exception as exc:  # network failures are reported, not hidden
            last_error = exc
            time.sleep(2**attempt)
    raise RuntimeError(f"Unable to retrieve {url}: {last_error}")


def ncbi_json(endpoint: str, params: dict[str, str]) -> dict:
    url = f"{EUTILS}/{endpoint}?{urllib.parse.urlencode(params)}"
    return json.loads(fetch(url).decode("utf-8"))


def read_dbids() -> list[str]:
    workbook = load_workbook(SOURCE, read_only=True, data_only=True)
    sheet = workbook.active
    header = [cell.value for cell in next(sheet.iter_rows())]
    expected = ["DBid", "CK-1", "CK-2", "CK-3", "LS-1", "LS-2", "LS-3", "MS-1", "MS-2", "MS-3"]
    if header != expected:
        raise ValueError(f"Unexpected table 6 header: {header}")
    dbids = [str(row[0].value) for row in sheet.iter_rows(min_row=2)]
    if len(dbids) != 100 or len(set(dbids)) != 100:
        raise ValueError("Expected exactly 100 unique DBids")
    return dbids


def blank_row(dbid: str, identifier_type: str) -> dict[str, str]:
    return {
        "DBid": dbid,
        "identifier_type": identifier_type,
        "resolved": "FALSE",
        "source_database": "",
        "organism": "",
        "gene_name": "",
        "protein_function_description": "",
        "KO_identifier": "",
        "EC_number": "",
        "evidence_source_URL_or_accession": "",
        "mapping_confidence": "low",
        "clearly_nitrogen_related": "FALSE",
        "Figure_7_category": "",
        "mapping_notes": "",
    }


def taxonomy_names(taxids: list[str]) -> dict[str, str]:
    if not taxids:
        return {}
    payload = ncbi_json("esummary.fcgi", {"db": "taxonomy", "id": ",".join(taxids), "retmode": "json"})
    result = payload.get("result", {})
    names: dict[str, str] = {}
    for uid in result.get("uids", []):
        record = result.get(uid, {})
        names[uid] = record.get("scientificname", "")
    return names


def resolve_gene_ids(rows: dict[str, dict[str, str]], gene_ids: list[str]) -> None:
    if not gene_ids:
        return
    payload = ncbi_json("esummary.fcgi", {"db": "gene", "id": ",".join(gene_ids), "retmode": "json"})
    result = payload.get("result", {})
    for gene_id in gene_ids:
        record = result.get(gene_id)
        if not record or record.get("error"):
            rows[gene_id]["mapping_notes"] = "No exact NCBI Gene record recovered."
            continue
        organism = record.get("organism", {})
        rows[gene_id].update(
            {
                "identifier_type": "bare numeric identifier; candidate NCBI Gene ID",
                "resolved": "FALSE",
                "source_database": "NCBI Gene (candidate match only)",
                "organism": organism.get("scientificname", ""),
                "gene_name": record.get("name", ""),
                "protein_function_description": record.get("description", ""),
                "evidence_source_URL_or_accession": f"https://www.ncbi.nlm.nih.gov/gene/{gene_id}",
                "mapping_confidence": "medium",
                "mapping_notes": "The bare integer matches this NCBI Gene record, but no database namespace is supplied; the identity cannot be confirmed uniquely and is not counted as resolved.",
            }
        )


def genpept_fields(uid: str) -> tuple[str, str, str]:
    url = f"{EUTILS}/efetch.fcgi?{urllib.parse.urlencode({'db': 'protein', 'id': uid, 'rettype': 'gb', 'retmode': 'text'})}"
    text = fetch(url).decode("utf-8", errors="replace")
    gene_match = re.search(r'/gene="([^"]+)"', text)
    product_match = re.search(r'/product="([^"]+(?:\n\s+[^"]+)*)"', text)
    product = re.sub(r"\s+", " ", product_match.group(1)).strip() if product_match else ""
    ecs = sorted(set(re.findall(r'/EC_number="([^"]+)"', text)))
    return (gene_match.group(1) if gene_match else "", product, ";".join(ecs))


def resolve_composite(rows: dict[str, dict[str, str]], dbid: str, taxon_name: str) -> str:
    taxid, locus = dbid.split(".", 1)
    term = f'"{locus}"[All Fields] AND txid{taxid}[Organism:exp]'
    search = ncbi_json("esearch.fcgi", {"db": "protein", "term": term, "retmode": "json", "retmax": "20"})
    ids = search.get("esearchresult", {}).get("idlist", [])
    if not ids:
        rows[dbid].update(
            {
                "source_database": "NCBI Taxonomy (prefix only)",
                "organism": taxon_name,
                "evidence_source_URL_or_accession": f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={taxid}",
                "mapping_notes": "Taxon prefix resolved, but no exact NCBI Protein record was recovered for the locus tag; full DBid remains unresolved.",
            }
        )
        return ""

    summaries = ncbi_json("esummary.fcgi", {"db": "protein", "id": ",".join(ids), "retmode": "json"}).get("result", {})
    matches = []
    for uid in ids:
        record = summaries.get(uid, {})
        if str(record.get("taxid", "")) == taxid and locus in record.get("extra", ""):
            matches.append((uid, record))
    if len(matches) != 1:
        rows[dbid].update(
            {
                "source_database": "NCBI Taxonomy (prefix only)",
                "organism": taxon_name,
                "evidence_source_URL_or_accession": f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={taxid}",
                "mapping_notes": f"Taxon prefix resolved, but the locus search produced {len(matches)} exact candidate protein records; no automatic choice was made.",
            }
        )
        return ""

    uid, record = matches[0]
    accession = record.get("accessionversion", "") or record.get("caption", "")
    gene, product, ec = genpept_fields(uid)
    description = product or re.sub(r"\s*\[[^\]]+\]\s*$", "", record.get("title", ""))
    rows[dbid].update(
        {
            "resolved": "TRUE",
            "source_database": "NCBI Protein; NCBI Taxonomy",
            "organism": record.get("organism", "") or taxon_name,
            "gene_name": gene,
            "protein_function_description": description,
            "EC_number": ec,
            "evidence_source_URL_or_accession": f"{accession}; https://www.ncbi.nlm.nih.gov/protein/{accession}",
            "mapping_confidence": "high",
            "mapping_notes": "Exact taxon ID and locus-tag match to one NCBI Protein record.",
        }
    )
    return accession


def kegg_ko_mapping(accessions: list[str]) -> dict[str, list[str]]:
    """Attempt direct KEGG conversion; absence of a result is not inferred."""
    mapping: dict[str, list[str]] = defaultdict(list)
    for start in range(0, len(accessions), 10):
        chunk = accessions[start : start + 10]
        query = "+".join(f"ncbi-proteinid:{x}" for x in chunk)
        url = f"https://rest.kegg.jp/conv/ko/{query}"
        text = fetch(url).decode("utf-8", errors="replace")
        for line in text.splitlines():
            if "\t" not in line:
                continue
            source, ko = line.split("\t", 1)
            accession = source.removeprefix("ncbi-proteinid:")
            mapping[accession].append(ko.removeprefix("ko:"))
    return mapping


def direct_nitrogen_classification(row: dict[str, str]) -> tuple[bool, str]:
    """Assign only categories stated unambiguously by recovered annotation."""
    text = f"{row['gene_name']} {row['protein_function_description']}".lower()
    if not text.strip():
        return False, ""

    nitrogen_terms = (
        "nitrogenase", "nitrate reductase", "nitrite reductase", "nitric oxide reductase",
        "nitrous oxide reductase", "nitrous-oxide reductase", "glutamine synthetase",
        "ammonia monooxygenase", "hydroxylamine",
        "nitrite oxidoreductase", "hydrazine synthase", "hydrazine dehydrogenase",
        "nitrogen fixation", "denitrification", "nitrification", "anammox",
    )
    related = any(term in text for term in nitrogen_terms)

    if "anammox" in text or "hydrazine synthase" in text or "hydrazine dehydrogenase" in text:
        return True, "Anammox"
    if "nitrification" in text or "ammonia monooxygenase" in text or "nitrite oxidoreductase" in text:
        return True, "Nitrification"
    if "nitrogen fixation" in text or "nitrogenase" in text:
        return True, "Nitrogen fixation"
    if (
        "denitrification" in text
        or "nitric oxide reductase" in text
        or "nitrous oxide reductase" in text
        or "nitrous-oxide reductase" in text
        or row["gene_name"].lower() == "nosz"
    ):
        return True, "Denitrification"
    if "glutamine synthetase" in text:
        return True, "Organic degradation and synthesis"
    if "respiratory nitrate reductase" in text or "respiratory nitrite reductase" in text:
        return True, "Dissimilatory nitrate reduction"
    if "dissimilatory nitrate reduction" in text or "dissimilatory nitrite reduction" in text:
        return True, "Dissimilatory nitrate reduction"
    if "assimilatory nitrate reduction" in text or "assimilatory nitrite reduction" in text:
        return True, "Assimilatory nitrate reduction"
    return related, ""


def write_mapping(rows: list[dict[str, str]]) -> None:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_audit(rows: list[dict[str, str]], network_notes: list[str]) -> None:
    total = len(rows)
    resolved = sum(row["resolved"] == "TRUE" for row in rows)
    with_function = sum(bool(row["gene_name"] or row["protein_function_description"]) for row in rows)
    resolved_with_function = sum(
        row["resolved"] == "TRUE" and bool(row["gene_name"] or row["protein_function_description"])
        for row in rows
    )
    candidate_function = with_function - resolved_with_function
    with_ko = sum(bool(row["KO_identifier"]) for row in rows)
    nitrogen = sum(row["clearly_nitrogen_related"] == "TRUE" for row in rows)
    category_rows: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        if row["Figure_7_category"]:
            category_rows[row["Figure_7_category"]].append(row["DBid"])

    lines = [
        "# Audit of nitrogen-related DBid recoverability",
        "",
        "## Scope and approach",
        "",
        "Supplementary `table 6.xlsx` was read without modification. Its single worksheet contains 100 unique DBid rows and nine abundance columns (`CK-1` through `MS-3`). This audit tests only whether each DBid can be mapped automatically and reproducibly to an authoritative public record. It does not infer annotations from identifier appearance or abundance patterns.",
        "",
        "Bare integer identifiers were tested against NCBI Gene, but matching records were retained as candidate-only because the supplement supplies no database namespace. Composite identifiers of the form `numeric prefix.locus tag` were interpreted syntactically as a candidate NCBI Taxonomy ID plus locus tag, then counted as resolved only when the locus tag and taxon matched one unique NCBI Protein record. Resolving the taxon prefix alone did not resolve the full DBid. Direct KEGG KO links were attempted for recovered NCBI protein accessions. Opaque 32-character hexadecimal strings were retained as unresolved because the supplementary material supplies no namespace or cross-reference table.",
        "",
        "Authoritative services used: [NCBI Entrez E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/), [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/), [NCBI Protein](https://www.ncbi.nlm.nih.gov/protein/), [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy), and the [KEGG REST API](https://www.kegg.jp/kegg/rest/keggapi.html). The per-row accessions and URLs are recorded in the mapping table.",
        "",
        "## Recovery summary",
        "",
        f"- Resolved full DBids: **{resolved}/{total} ({resolved / total:.1%})**",
        f"- Rows with a recovered gene name and/or protein/function description: **{with_function}/{total} ({with_function / total:.1%})**; this comprises **{resolved_with_function}** reliable full-identifier mappings and **{candidate_function}** candidate-only bare-integer matches",
        f"- Rows with a direct KO mapping: **{with_ko}/{total} ({with_ko / total:.1%})**",
        f"- Rows clearly nitrogen-related from the recovered annotation: **{nitrogen}/{total} ({nitrogen / total:.1%})**",
        f"- Unresolved full DBids: **{total - resolved}/{total} ({(total - resolved) / total:.1%})**",
        "",
        "Identifier classes:",
        "",
        "| Identifier class | Rows | Recoverability |",
        "|---|---:|---|",
    ]
    final_class_counts = Counter(row["identifier_type"] for row in rows)
    for cls, count in sorted(final_class_counts.items()):
        class_resolved = sum(row["identifier_type"] == cls and row["resolved"] == "TRUE" for row in rows)
        lines.append(f"| {cls} | {count} | {class_resolved}/{count} full identifiers resolved |")

    lines.extend([
        "",
        "## Comparison with Figure 7 categories",
        "",
        "The category comparison below is intentionally strict. A category is recorded only when the recovered gene/protein annotation explicitly supports it. Generic nitrate- or nitrite-reductase wording is nitrogen-related but is not sufficient by itself to distinguish assimilatory reduction, dissimilatory reduction, or denitrification. `Organic degradation and synthesis` is assigned only to explicit glutamine-synthetase annotations; the broader contents of that paper-defined category cannot be reconstructed. `Auxiliary` is not assigned automatically because the paper does not publish a reproducible DBid-to-category rule for it.",
        "",
        "| Figure 7 category | Directly supported DBids | Count |",
        "|---|---|---:|",
    ])
    for category in FIGURE7_CATEGORIES:
        identifiers = category_rows.get(category, [])
        lines.append(f"| {category} | {'; '.join(identifiers) if identifiers else 'None'} | {len(identifiers)} |")

    lines.extend([
        "",
        "## Reproducibility assessment",
        "",
        "The recoverable records provide evidence for individual gene or protein annotations, but they do not reconstruct the paper's complete mapping of all 100 DBids to the nine Figure 7 categories. In particular, opaque hashes cannot be linked to genes without the authors' original annotation table or gene catalogue, and broad functional categories cannot be recreated from abundance values alone. Consequently, Figure 7 cannot be reproduced faithfully from `table 6.xlsx` alone even where some individual identifiers resolve.",
        "",
        "The mapping output preserves unresolved rows rather than substituting similarity-based guesses. A blank KO, EC number, or Figure 7 category means that no direct mapping was recovered; it does not demonstrate that the underlying sequence lacks that annotation.",
    ])
    if network_notes:
        lines.extend(["", "## Retrieval limitations", "", *[f"- {note}" for note in network_notes]])
    lines.extend([
        "",
        "## Outputs",
        "",
        "- Mapping table: `data/processed/nitrogen_dbid_mapping.tsv`",
        "- Reproducible resolver: `scripts/resolve_nitrogen_dbids.py`",
        "",
        "The resolver requires the original `table 6.xlsx` at `references/lv_2026/table 6.xlsx`, internet access to the named services and Python with `openpyxl` (declared in `workflow/environment.yml`). Run it from the project root with the `alfalfa-snakemake` Python interpreter. It refreshes the mapping TSV and this audit; live database records can change, so the deposited mapping remains the reference for the published findings. Ordinary analytical reproduction uses the validated processed inputs and does not require repeating these remote queries.",
        "",
    ])
    AUDIT.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    dbids = read_dbids()
    rows: dict[str, dict[str, str]] = {}
    gene_ids: list[str] = []
    composite: list[tuple[str, str, str]] = []

    for dbid in dbids:
        if HEX_RE.fullmatch(dbid):
            identifier_type = "opaque 32-character hexadecimal identifier"
            rows[dbid] = blank_row(dbid, identifier_type)
            rows[dbid]["mapping_notes"] = "No identifier namespace or public cross-reference is supplied; not guessed."
        elif INTEGER_RE.fullmatch(dbid):
            identifier_type = "candidate NCBI Gene ID"
            rows[dbid] = blank_row(dbid, identifier_type)
            gene_ids.append(dbid)
        else:
            match = COMPOSITE_RE.fullmatch(dbid)
            if not match:
                identifier_type = "unclassified identifier"
                rows[dbid] = blank_row(dbid, identifier_type)
                rows[dbid]["mapping_notes"] = "Identifier syntax was not recognized; not guessed."
            else:
                identifier_type = "NCBI Taxonomy ID plus locus tag"
                rows[dbid] = blank_row(dbid, identifier_type)
                composite.append((dbid, match.group(1), match.group(2)))
    resolve_gene_ids(rows, gene_ids)
    taxids = sorted({taxid for _, taxid, _ in composite}, key=int)
    taxa = taxonomy_names(taxids)
    accession_to_dbid: dict[str, str] = {}
    for index, (dbid, taxid, _locus) in enumerate(composite, 1):
        print(f"Resolving locus {index}/{len(composite)}: {dbid}", flush=True)
        accession = resolve_composite(rows, dbid, taxa.get(taxid, ""))
        if accession:
            accession_to_dbid[accession] = dbid

    network_notes: list[str] = []
    try:
        ko_map = kegg_ko_mapping(list(accession_to_dbid))
        for accession, kos in ko_map.items():
            dbid = accession_to_dbid.get(accession)
            if dbid:
                rows[dbid]["KO_identifier"] = ";".join(sorted(set(kos)))
    except Exception as exc:
        network_notes.append(f"Direct KEGG KO conversion could not be completed: {exc}")

    ordered_rows = [rows[dbid] for dbid in dbids]
    for row in ordered_rows:
        related, category = direct_nitrogen_classification(row)
        row["clearly_nitrogen_related"] = "TRUE" if related else "FALSE"
        row["Figure_7_category"] = category

    write_mapping(ordered_rows)
    write_audit(ordered_rows, network_notes)
    print(f"Wrote {OUTPUT}")
    print(f"Wrote {AUDIT}")


if __name__ == "__main__":
    main()
