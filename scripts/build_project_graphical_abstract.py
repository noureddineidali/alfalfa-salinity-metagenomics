#!/usr/bin/env python3
"""Compose an editable, self-contained graphical abstract from validated figures.

The three source SVG plots remain native vector elements in the output. The
validated functional PCoA is embedded unchanged from its published PNG export.
This script performs no ordination, filtering, statistics, or data conversion.
"""

import base64
import hashlib
from pathlib import Path
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/figures/project_graphical_abstract.svg"
SOURCES = {
    "genus": ROOT / "results/figures/report1-presentation/genus_bray_curtis_pcoa_refined.svg",
    "functional": ROOT / "results/figures/normalized_functional_pathway_pcoa.png",
    "read_fate": ROOT / "results/figures/raw-read-report/read_fate_ck1_ms1.svg",
    "comparison": ROOT / "results/figures/raw-read-report/published_vs_reconstructed_ck1_ms1.svg",
}


def native_svg(path: Path, x: int, y: int, width: int, height: int) -> str:
    root = ET.parse(path).getroot()
    view_box = root.attrib["viewBox"]
    contents = "".join(ET.tostring(child, encoding="unicode") for child in root)
    return (
        f'<svg x="{x}" y="{y}" width="{width}" height="{height}" '
        f'viewBox="{view_box}" preserveAspectRatio="xMidYMid meet">'
        f"{contents}</svg>"
    )


functional_uri = "data:image/png;base64," + base64.b64encode(SOURCES["functional"].read_bytes()).decode("ascii")
source_manifest = "\n".join(
    f"    {path.relative_to(ROOT)}  SHA-256 {hashlib.sha256(path.read_bytes()).hexdigest()}"
    for path in SOURCES.values()
)

svg = f'''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"
     width="2200" height="1280" viewBox="0 0 2200 1280" role="img"
     aria-labelledby="abstract-title abstract-description">
  <title id="abstract-title">Two evidence streams for the alfalfa root-zone microbiome</title>
  <desc id="abstract-description">A stylized alfalfa plant and salinity gradient introduce two complementary reanalyses. Report 1 shows the actual nine-sample genus and normalized functional ordinations for CK, LS and MS. Report 2 shows the actual CK-1 and MS-1 read-fate and published-versus-reconstructed taxonomic comparison figures. Findings are descriptive within the limits of replication and database coverage.</desc>
  <!-- Rebuilt from the previous editable project abstract; plot coordinates and counts are not redrawn or simulated.
  Validated source figures:
{source_manifest}
  -->
  <defs>
    <linearGradient id="hero-wash" x1="0" y1="0" x2="1" y2="0">
      <stop offset="0" stop-color="#e9f4ef"/><stop offset=".58" stop-color="#edf5f4"/><stop offset="1" stop-color="#f6eee4"/>
    </linearGradient>
    <linearGradient id="soil-wash" x1="0" y1="0" x2="1" y2="0">
      <stop offset="0" stop-color="#b8aea0"/><stop offset="1" stop-color="#c7b59f"/>
    </linearGradient>
    <style>
      .ink{{fill:#263238}}.muted{{fill:#5f6b72}}.teal{{fill:#007e72}}.coral{{fill:#be684d}}
      .eyebrow{{font:700 22px Inter,Arial,sans-serif;letter-spacing:2.5px}}
      .hero-title{{font:750 53px Inter,Arial,sans-serif;letter-spacing:-1.1px}}
      .hero-sub{{font:400 26px Inter,Arial,sans-serif}}
      .card-title{{font:700 32px Inter,Arial,sans-serif;letter-spacing:-.4px}}
      .card-sub{{font:400 21px Inter,Arial,sans-serif}}
      .figure-label{{font:700 22px Inter,Arial,sans-serif}}
      .finding{{font:600 23px Inter,Arial,sans-serif}}
      .small{{font:400 19px Inter,Arial,sans-serif}}
      .footer-main{{font:650 24px Inter,Arial,sans-serif}}
    </style>
  </defs>

  <rect width="2200" height="1280" fill="#f6f9fa"/>
  <rect width="2200" height="15" fill="#263e4e"/>
  <rect x="0" y="15" width="2200" height="268" fill="url(#hero-wash)"/>

  <!-- Botanical context: vector leaflets, stem, root branches and soil. No nodules or stress phenotype is implied. -->
  <path d="M0 205 C350 192 690 219 1070 201 S1790 211 2200 195 L2200 283 L0 283Z" fill="url(#soil-wash)" opacity=".82"/>
  <path d="M0 204 C350 191 690 218 1070 200 S1790 210 2200 194" fill="none" stroke="#8d9d88" stroke-width="3"/>
  <g opacity=".4" fill="#857b70">
    <circle cx="84" cy="245" r="4"/><circle cx="180" cy="226" r="3"/><circle cx="374" cy="251" r="4"/>
    <circle cx="470" cy="224" r="3"/><circle cx="700" cy="249" r="3"/><circle cx="920" cy="233" r="4"/>
    <circle cx="1308" cy="253" r="4"/><circle cx="1580" cy="222" r="3"/><circle cx="1880" cy="246" r="4"/>
  </g>
  <g fill="none" stroke="#6e665a" stroke-width="3" stroke-linecap="round" opacity=".8">
    <path d="M285 199 C277 222 268 238 255 263"/><path d="M285 199 C299 219 309 235 324 261"/>
    <path d="M281 214 C265 222 245 226 229 241"/><path d="M288 218 C305 224 324 225 344 239"/>
    <path d="M272 232 C269 244 264 252 260 272"/><path d="M309 237 C319 246 326 254 331 272"/>
  </g>
  <g fill="none" stroke="#317f5e" stroke-width="6" stroke-linecap="round">
    <path d="M284 200 C284 169 286 132 279 91"/>
    <path d="M284 154 C263 144 240 135 221 119"/><path d="M284 132 C306 117 329 107 343 92"/>
    <path d="M281 112 C265 96 252 81 242 64"/>
  </g>
  <g fill="#3f9b6b" stroke="#2f865f" stroke-width="2">
    <path d="M219 119 C194 100 178 88 172 72 C191 72 212 84 226 108Z"/>
    <path d="M218 116 C206 90 207 72 215 57 C230 68 235 91 225 111Z"/>
    <path d="M222 116 C242 98 258 93 273 96 C265 113 246 123 225 123Z"/>
    <path d="M343 92 C357 67 373 60 390 59 C387 79 370 94 346 101Z"/>
    <path d="M342 92 C337 68 341 51 352 38 C365 54 364 74 348 96Z"/>
    <path d="M343 93 C365 94 382 105 389 121 C368 120 350 110 340 99Z"/>
    <path d="M245 67 C231 46 229 29 236 17 C251 29 257 48 250 69Z"/>
    <path d="M245 67 C254 43 268 34 285 31 C278 50 262 63 248 72Z"/>
    <path d="M245 67 C224 61 208 62 195 70 C211 80 230 78 247 71Z"/>
  </g>
  <g fill="none" stroke="#b58d5e" stroke-width="3" opacity=".72">
    <path d="M421 233 l9 -9 9 9 -9 9Z"/><path d="M487 255 l6 -6 6 6 -6 6Z"/>
    <path d="M382 261 l5 -5 5 5 -5 5Z"/>
  </g>

  <text x="545" y="73" class="eyebrow teal">PRJNA1399466  ·  TWO COMPLEMENTARY REANALYSES</text>
  <text x="545" y="135" class="hero-title ink">Alfalfa root-zone metagenomics under salinity</text>
  <text x="547" y="178" class="hero-sub muted">Published profiles reveal community patterns; deposited reads test taxonomic reproducibility.</text>
  <g font-family="Inter,Arial,sans-serif" font-size="20" font-weight="700">
    <rect x="545" y="211" width="208" height="40" rx="20" fill="#e9f2f8"/><circle cx="568" cy="231" r="7" fill="#0072b2"/><text x="585" y="238" fill="#165676">CK · 0 g kg⁻¹</text>
    <rect x="770" y="211" width="233" height="40" rx="20" fill="#fff3d9"/><circle cx="793" cy="231" r="7" fill="#e69f00"/><text x="810" y="238" fill="#866100">LS · 3.0 g kg⁻¹</text>
    <rect x="1020" y="211" width="238" height="40" rx="20" fill="#e6f3eb"/><circle cx="1043" cy="231" r="7" fill="#009e73"/><text x="1060" y="238" fill="#1a7055">MS · 4.5 g kg⁻¹</text>
  </g>

  <!-- Two analytical branches; the earlier abstract's paired-card composition is retained. -->
  <rect x="62" y="314" width="1017" height="824" rx="20" fill="#fff" stroke="#d8e3e7" stroke-width="2"/>
  <rect x="1121" y="314" width="1017" height="824" rx="20" fill="#fff" stroke="#d8e3e7" stroke-width="2"/>
  <rect x="62" y="314" width="1017" height="9" rx="4" fill="#008e70"/>
  <rect x="1121" y="314" width="1017" height="9" rx="4" fill="#c96a48"/>

  <text x="94" y="365" class="eyebrow teal">01  /  PUBLISHED MATRICES</text>
  <text x="94" y="409" class="card-title ink">Nine-sample ecological reanalysis</text>
  <text x="94" y="441" class="card-sub muted">CK, LS and MS · three sequenced samples per treatment</text>
  <rect x="92" y="466" width="468" height="416" rx="13" fill="#fbfcfd" stroke="#e3ebee"/>
  <rect x="581" y="466" width="468" height="416" rx="13" fill="#fbfcfd" stroke="#e3ebee"/>
  {native_svg(SOURCES["genus"], 99, 490, 456, 345)}
  <image x="591" y="493" width="448" height="345" href="{functional_uri}" preserveAspectRatio="xMidYMid meet"/>
  <text x="108" y="863" class="figure-label ink">Genus Bray–Curtis PCoA</text>
  <text x="596" y="863" class="figure-label ink">Normalized pathway PCoA</text>
  <text x="94" y="922" class="finding ink">Community structure: R² = 0.775  ·  alpha diversity overlaps</text>
  <text x="94" y="965" class="finding ink">Relative 30-pathway structure: R² = 0.822</text>
  <text x="94" y="1007" class="small muted">Raw functional PCoA1 tracks pathway total (r = −0.9988).</text>
  <text x="94" y="1043" class="small muted">CK-3 differs in composition; measured soil values remain CK-like.</text>
  <rect x="92" y="1070" width="957" height="42" rx="8" fill="#ecf6f1"/>
  <text x="110" y="1098" class="small" fill="#326754">Genomic potential is not activity  ·  nitrogen Figure 7 lacks a full identifier crosswalk</text>

  <text x="1153" y="365" class="eyebrow coral">02  /  DEPOSITED SRA READS</text>
  <text x="1153" y="409" class="card-title ink">Independent CK-1 / MS-1 reconstruction</text>
  <text x="1153" y="441" class="card-sub muted">SRA → QC → host depletion → Kraken2 / Bracken</text>
  <rect x="1151" y="466" width="957" height="295" rx="13" fill="#fbfcfd" stroke="#e3ebee"/>
  {native_svg(SOURCES["read_fate"], 1164, 474, 930, 276)}
  <rect x="1151" y="778" width="957" height="295" rx="13" fill="#fbfcfd" stroke="#e3ebee"/>
  {native_svg(SOURCES["comparison"], 1164, 787, 930, 276)}
  <text x="1153" y="1105" class="small muted">Both runs validated; exact-label disagreement does not establish taxon absence.</text>

  <rect x="62" y="1167" width="2076" height="82" rx="16" fill="#263e4e"/>
  <text x="92" y="1202" class="eyebrow" fill="#a9d8d4">READ TOGETHER</text>
  <text x="92" y="1232" class="footer-main" fill="#fff">Ecological patterns and reconstruction choices both matter. Replication, reference coverage and incomplete annotations limit interpretation.</text>
</svg>
'''

OUT.write_text(svg, encoding="utf-8")
ET.parse(OUT)
print(f"Wrote {OUT.relative_to(ROOT)} from four validated source figures")
