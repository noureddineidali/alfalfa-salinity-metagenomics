"""Assemble and check the small, public GitHub Pages artifact."""

from __future__ import annotations

import argparse
import html
import shutil
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit


ROOT = Path(__file__).resolve().parents[1]
PUBLIC_FILES = (
    "report/alfalfa_salinity_metagenome_reanalysis.html",
    "report/alfalfa_ck1_raw_read_reconstruction.html",
    "results/figures/project_graphical_abstract.svg",
    "docs/dataset_reconnaissance.md",
    "docs/supplementary_data_audit.md",
    "docs/processed_data_validation.md",
    "docs/nitrogen_mapping_audit.md",
    "docs/raw_read_provenance.tsv",
    "workflow/README.md",
    "results/workflow/raw_read_scientific_overview.svg",
    "results/workflow/raw_read_scientific_overview.dot",
    "results/workflow/ms1_raw_read_workflow.svg",
    "results/workflow/ms1_raw_read_workflow.dot",
    "results/workflow/raw_read_rule_graph.svg",
    "results/workflow/raw_read_rule_graph.dot",
    "data/metadata/sample_manifest.tsv",
    "data/processed/genus_abundance.tsv",
)
SUPPORT_DIRECTORIES = (
    "analysis/raw_read_comparison",
    "results/validation/raw-read-workflow/MS-1",
    "results/benchmarks/raw-read-workflow/MS-1",
)


class Links(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.local: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        for key, value in attrs:
            if key in {"href", "src"} and value:
                parsed = urlsplit(value)
                if not parsed.scheme and not parsed.netloc and parsed.path:
                    self.local.append(unquote(parsed.path))


def copy_file(destination: Path, relative: str) -> None:
    source = ROOT / relative
    if not source.is_file():
        raise FileNotFoundError(source)
    target = destination / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def directory_index(title: str, files: list[Path]) -> str:
    items = "\n".join(
        f'<li><a href="{html.escape(path.name, quote=True)}">'
        f"{html.escape(path.name)}</a></li>"
        for path in files
    )
    return (
        "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">"
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        f"<title>{html.escape(title)}</title>"
        "<style>body{font:16px/1.6 system-ui,sans-serif;max-width:760px;"
        "margin:4rem auto;padding:0 1.5rem;color:#17313b}a{color:#205f83}"
        "li{margin:.5rem 0}</style></head><body>"
        f"<h1>{html.escape(title)}</h1><ul>{items}</ul></body></html>\n"
    )


def build(destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / "site/index.html", destination / "index.html")
    (destination / ".nojekyll").touch()
    for relative in PUBLIC_FILES:
        copy_file(destination, relative)
    for relative in SUPPORT_DIRECTORIES:
        source = ROOT / relative
        files = sorted(path for path in source.iterdir() if path.is_file())
        if not files:
            raise ValueError(f"No public support files in {source}")
        for path in files:
            copy_file(destination, f"{relative}/{path.name}")
        target = destination / relative / "index.html"
        target.write_text(directory_index(relative.split("/")[-1], files), encoding="utf-8")


def validate(destination: Path) -> None:
    missing: list[str] = []
    for page in destination.rglob("*.html"):
        parser = Links()
        parser.feed(page.read_text(encoding="utf-8"))
        for link in parser.local:
            target = destination / link.lstrip("/") if link.startswith("/") else page.parent / link
            if target.is_dir():
                target /= "index.html"
            if not target.exists():
                missing.append(f"{page.relative_to(destination)}: {link}")
    if missing:
        raise ValueError("Unresolved site links/assets:\n" + "\n".join(missing))
    count = sum(1 for path in destination.rglob("*") if path.is_file())
    print(f"Pages artifact validated: {count} files; all local HTML links and assets resolve")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    build(args.destination)
    validate(args.destination)
