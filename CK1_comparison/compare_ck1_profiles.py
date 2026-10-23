#!/usr/bin/env python3
"""Reproduce exact-label comparison of Lv et al. table-4 CK-1 vs PlusPF-8 Bracken.

Requirements: Python 3, numpy, scipy, matplotlib.
Usage: python compare_ck1_profiles.py --published ck1_published_table4_extracted.csv \
        --bracken CK1_pluspf8_c0_genus.bracken --outdir output
This script does not synonym-normalize taxa or infer taxon equivalence.
"""
import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr, spearmanr
from scipy.spatial.distance import braycurtis, jensenshannon


def write_csv(path, fields, rows):
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def fmt(v, digits=2):
    return f'{v:.{digits}f}'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--published', type=Path, required=True)
    ap.add_argument('--bracken', type=Path, required=True)
    ap.add_argument('--outdir', type=Path, required=True)
    args = ap.parse_args()
    out = args.outdir
    out.mkdir(exist_ok=True, parents=True)
    with args.published.open(newline='', encoding='utf-8') as f:
        p_rows = list(csv.DictReader(f))
    published = {r['genus_label'].strip(): float(r['published_original_abundance']) for r in p_rows}
    assert len(published) == len(p_rows) == 15805
    with args.bracken.open(newline='', encoding='utf-8') as f:
        b_rows = list(csv.DictReader(f, delimiter='\t'))
    bracken_rows = {r['name'].strip(): r for r in b_rows}
    bracken = {k: int(v['new_est_reads']) for k,v in bracken_rows.items()}
    assert len(bracken)==len(b_rows)==1897
    assert all(r['taxonomy_lvl']=='G' for r in b_rows)
    Ptotal = sum(published.values())
    Btotal = sum(bracken.values())
    P={k:v/Ptotal for k,v in published.items()}
    B={k:v/Btotal for k,v in bracken.items()}
    names=set(published)|set(bracken)
    both=set(published)&set(bracken)
    both_pos={x for x in both if published[x]>0 and bracken[x]>0}
    p_sort=sorted(P,key=lambda x:(-P[x],x))
    b_sort=sorted(B,key=lambda x:(-B[x],x))
    ranks_p={x:i+1 for i,x in enumerate(p_sort)}
    ranks_b={x:i+1 for i,x in enumerate(b_sort)}
    intersect={k: len(set(p_sort[:k]) & set(b_sort[:k])) for k in [10,20,50,100,500]}
    ks=sorted(both_pos)
    x=np.array([P[z] for z in ks]); y=np.array([B[z] for z in ks])
    union_names=sorted(names)
    pu=np.array([P.get(z,0.0) for z in union_names]); bu=np.array([B.get(z,0.0) for z in union_names]);
    spearman=spearmanr(x,y)
    pearson_logs=pearsonr(np.log10(x),np.log10(y))
    bc=float(braycurtis(pu,bu))
    jsd=float(jensenshannon(pu,bu,base=2)**2)
    metrics={
        'published_genus_labels':len(published),
        'published_nonzero_labels': sum(1 for v in published.values() if v>0),
        'bracken_genus_labels':len(bracken),
        'exact_shared_labels':len(both),
        'shared_labels_positive_in_both':len(both_pos),
        'published_labels_positive_not_in_bracken':sum(v>0 and k not in bracken for k,v in published.items()),
        'bracken_labels_absent_from_published':len(set(bracken)-set(published)),
        'published_abundance_share_of_shared_labels':float(sum(P[k] for k in both)),
        'bracken_abundance_share_of_shared_labels':float(sum(B[k] for k in both)),
        'spearman_shared_positive':float(spearman.statistic),
        'pearson_log10_shared_positive':float(pearson_logs.statistic),
        'bray_curtis_union_exact_labels':bc,
        'jensen_shannon_divergence_base2_union_exact_labels':jsd,
        'topk_intersections':{str(k):v for k,v in intersect.items()},
        'published_original_abundance_total_unknown_units':Ptotal,
        'bracken_sum_new_est_reads':Btotal,
        'bracken_sum_kraken_assigned_reads':sum(int(r['kraken_assigned_reads']) for r in b_rows),
        'bracken_sum_added_reads':sum(int(r['added_reads']) for r in b_rows),
        'bracken_sum_fraction_total_reads':sum(float(r['fraction_total_reads']) for r in b_rows),
        'kraken2_classified_pairs_from_session_log':1833747,
        'kraken2_total_pairs_from_session_log':33808912,
        'bracken_log_genus_level_reads_before_rowwise_integer_rounding':1830146,
        'bracken_log_vs_output_count_difference':1830146-Btotal,
    }
    with (out/'comparison_metrics.json').open('w') as f:
        json.dump(metrics,f,indent=2)

    all_rows=[]
    for name in union_names:
        pr=P.get(name)
        br=B.get(name)
        all_rows.append({
            'genus_label_exact':name,
            'match_status':'shared_label' if name in both else ('published_only' if name in P else 'bracken_only'),
            'published_original_abundance':published.get(name,''),
            'published_percent':100*pr if pr is not None else '',
            'bracken_new_est_reads':bracken.get(name,''),
            'bracken_percent_renormalized':100*br if br is not None else '',
            'bracken_reported_fraction_total_reads':bracken_rows[name]['fraction_total_reads'] if name in B else '',
            'bracken_taxid':bracken_rows[name]['taxonomy_id'] if name in B else '',
            'published_rank':ranks_p.get(name,''),
            'bracken_rank':ranks_b.get(name,''),
        })
    write_csv(out/'all_genus_exact_label_comparison.csv',list(all_rows[0]),all_rows)
    shared_rows=[r for r in all_rows if r['match_status']=='shared_label' and r['published_original_abundance']>0]
    write_csv(out/'shared_positive_genera.csv',list(all_rows[0]),shared_rows)

    toprows=[]
    for profile, ranking in [('published',p_sort),('bracken',b_sort)]:
        for i,nm in enumerate(ranking[:20],start=1):
            top=next(z for z in all_rows if z['genus_label_exact']==nm)
            toprows.append({'source_profile':profile,'rank_in_source':i,**top})
    write_csv(out/'top20_each_profile.csv',list(toprows[0]),toprows)

    # Plot styles: plain, readable, no inference about taxonomy synonyms.
    plt.rcParams.update({'font.size':10,'axes.titlesize':13,'axes.labelsize':10,'figure.dpi':170,'savefig.dpi':220})
    for label,rank,pr,filename in [('Published CK-1',p_sort,P,'fig1_published_top20.png'),
                                    ('Bracken PlusPF-8 CK-1',b_sort,B,'fig2_bracken_top20.png')]:
        selected=rank[:20][::-1]
        fig,ax=plt.subplots(figsize=(10,7.6))
        ax.barh(range(len(selected)),[pr[k]*100 for k in selected]);ax.set_yticks(range(len(selected)), selected)
        ax.set_xlabel('Relative abundance within source genus profile (%)')
        ax.set_title('Top 20 genus labels — '+label,loc='left',pad=12)
        ax.grid(axis='x',alpha=.17)
        ax.set_axisbelow(True)
        ax.spines[['top','right']].set_visible(False)
        fig.tight_layout()
        fig.savefig(out/filename,bbox_inches='tight')
        plt.close(fig)

    # Strong shared labels selected without selecting for differences (highest sum of source percentages).
    select=sorted(both_pos,key=lambda k:-(P[k]+B[k]))[:18]
    select=select[::-1]
    yy=np.arange(len(select))
    fig,ax=plt.subplots(figsize=(11,7.4))
    h=.38
    ax.barh(yy-h/2,[P[n]*100 for n in select],height=h,label='Lv et al. published CK-1')
    ax.barh(yy+h/2,[B[n]*100 for n in select],height=h,label='Independent Bracken CK-1')
    ax.set_yticks(yy, select)
    ax.set_xlabel('Within-profile relative abundance (%)')
    ax.set_title('Exact shared genera with the highest combined abundance',loc='left',pad=12)
    ax.legend(frameon=False,loc='lower right')
    ax.grid(axis='x',alpha=.17)
    ax.set_axisbelow(True);ax.spines[['top','right']].set_visible(False)
    fig.tight_layout();fig.savefig(out/'fig3_shared_genera_comparison.png',bbox_inches='tight');plt.close(fig)

    fig,ax=plt.subplots(figsize=(7.9,6.6))
    ax.scatter(x*100,y*100,s=12,alpha=.45,rasterized=True)
    lows=min(x.min(),y.min())*100*.7; highs=max(x.max(),y.max())*100*1.5
    ax.plot([lows,highs],[lows,highs],linestyle='--',linewidth=1,label='Equal abundance')
    ax.set_xscale('log');ax.set_yscale('log');ax.set_xlim(lows,highs);ax.set_ylim(lows,highs)
    ax.set_xlabel('Published CK-1 relative abundance (%) — log scale')
    ax.set_ylabel('Bracken CK-1 relative abundance (%) — log scale')
    ax.set_title(f'1,263 exact shared positive-abundance genus labels\nSpearman ρ = {spearman.statistic:.3f}',loc='left',pad=11)
    for tax in ['Bradyrhizobium','Streptomyces','Nocardioides','Mesorhizobium','Mycobacterium']:
        if tax in both_pos:
            ax.annotate(tax,(P[tax]*100,B[tax]*100),xytext=(3,4),textcoords='offset points',fontsize=8)
    ax.grid(True,which='major',alpha=.16);ax.legend(frameon=False)
    ax.spines[['top','right']].set_visible(False)
    fig.tight_layout();fig.savefig(out/'fig4_shared_positive_scatter.png',bbox_inches='tight');plt.close(fig)

    fig,ax=plt.subplots(figsize=(8,4.8))
    xr=np.array([10,20,50,100,200,500])
    intersections=[len(set(p_sort[:int(k)]) & set(b_sort[:int(k)])) for k in xr]
    ax.plot(xr,intersections,marker='o',label='Observed exact-name intersection')
    ax.set_xlabel('Top K genera considered from each profile')
    ax.set_ylabel('Number of exact labels shared within both top-K lists')
    ax.set_title('Top-K concordance across independently normalized profiles',loc='left',pad=12)
    ax.grid(True,alpha=.17);ax.spines[['top','right']].set_visible(False)
    fig.tight_layout();fig.savefig(out/'fig5_topk_overlap.png',bbox_inches='tight');plt.close(fig)

    # Standalone report for insertion into report 2; unlike a treatment comparison, this is one biological sample.
    top_b=sorted(bracken,key=lambda z:-B[z])[:10]
    top_p=sorted(published,key=lambda z:-P[z])[:10]
    def compare_row(n):
        pv=(f'{P[n]*100:.3f}%' if n in P else 'not listed')
        bv=(f'{B[n]*100:.3f}%' if n in B else 'not listed')
        return f'| {n} | {pv} | {bv} |'
    txt=f'''# CK-1 genus-profile comparison: Lv et al. (2026) versus Kraken2–Bracken PlusPF-8

**Study:** Lv et al. (2026), *Frontiers in Plant Science*, DOI: 10.3389/fpls.2026.1753229<br>
**Sample:** CK-1 (`SRR36739085`), the same biological sample in both datasets<br>
**Status:** Completed exploratory single-sample comparison; not a treatment-level replication.

## Sources and method

1. **Published:** `table 4.xlsx`, `Sheet1`, header `genus`, first numeric column (`1`) confirmed to represent CK-1 against the published `data sheet 2.csv` percentages. Contains {len(published):,} named genus-level labels. The published numerical unit is undocumented; original values were preserved and normalized by the CK-1 column sum ({Ptotal:,.6f} original units).
2. **Independent:** `CK1_pluspf8_c0_genus.bracken`, Bioconda Bracken package 3.1 (whose bundled executable reports v3.0.1; `-l G -r 150 -t 10`), using Kraken2 2.17.2 and the MD5-verified 26 June 2026 PlusPF-8 reference. Normalization for quantitative comparisons uses `new_est_reads / sum(new_est_reads)` ({Btotal:,} row-summed estimated fragments). Bracken's reported `fraction_total_reads` is also retained unchanged in the full CSV; the two representations differ slightly because abundance is calculated as floating point and exported as integer counts.
3. **Matching:** Exact, case-sensitive genus-label equality. No inference of taxonomic synonymy or equivalence for names from different taxonomic frameworks. Absence of a label in a source is marked as *not listed*, **not** as a confirmed biological zero.
4. **Statistics:** Descriptive genus-level comparison within **one sample**: Spearman rank correlation over labels with positive abundance in both sources, log10-scale Pearson correlation, top-K exact-label intersections, and Bray–Curtis/Jensen–Shannon between normalized vectors on the exact-label union. The distances are database-/nomenclature-dependent *profile differences*, not ecological between-sample diversity or treatment effect.

Twelve label cells in the published workbook were stored as Excel date-formatted numbers. The extracted input preserves their underlying serial values because their intended original genus strings cannot be recovered. None matches a Bracken genus; substituting the ISO date strings used in the processed repository table leaves all reported numerical metrics unchanged.

## Quantitative findings

| Measure | Value |
|---|---:|
| Labels in published genus matrix | {len(published):,} |
| Positive published CK-1 genus labels | {metrics['published_nonzero_labels']:,} |
| Genera in Bracken output | {len(bracken):,} |
| Labels exactly shared across both tables | {len(both):,} |
| Shared labels with positive abundance in both | {len(both_pos):,} |
| Bracken genera absent as exact labels from published table | {metrics['bracken_labels_absent_from_published']:,} |
| Positive published CK-1 labels absent as exact labels from Bracken | {metrics['published_labels_positive_not_in_bracken']:,} |
| Bracken abundance represented by exactly shared labels | {100*metrics['bracken_abundance_share_of_shared_labels']:.2f}% |
| Published abundance represented by exactly shared labels | {100*metrics['published_abundance_share_of_shared_labels']:.2f}% |
| Spearman ρ (shared positive, n={len(both_pos):,}) | {spearman.statistic:.3f} |
| Pearson r of log10 proportional abundances (same subset) | {pearson_logs.statistic:.3f} |
| Intersection of the two top-10 lists | {intersect[10]}/10 |
| Intersection of the two top-20 lists | {intersect[20]}/20 |
| Bray–Curtis on exact-label union | {bc:.3f} |
| Jensen–Shannon divergence (base 2) on exact-label union | {jsd:.3f} |

### Published CK-1 dominant labels

| Published rank | Genus label | Published CK-1 (%) | Bracken CK-1 (%) |
|---:|---|---:|---|
'''
    for i,k in enumerate(top_p,1):
        txt+=f'| {i} | {k} | {P[k]*100:.3f} | '+(f'{B[k]*100:.3f}' if k in B else 'not listed')+' |\n'
    txt+='''
### Reconstructed CK-1 dominant genera

| Bracken rank | Genus | Bracken CK-1 (%) | Published CK-1 (%) |
|---:|---|---:|---|
'''
    for i,k in enumerate(top_b,1):
        txt+=f'| {i} | {k} | {B[k]*100:.3f} | '+(f'{P[k]*100:.3f}' if k in P else 'not listed')+' |\n'
    txt+=f'''
## Interpretation

The independent reconstruction recovers a coherent set of soil-associated genera, including *Bradyrhizobium*, *Mesorhizobium*, *Streptomyces*, *Nocardioides*, *Sphingomonas* and *Arthrobacter*. Of {len(both):,} exactly shared genus labels, {len(both_pos):,} have positive abundance in both sources, and those shared positive labels show a moderately strong positive abundance-rank association (Spearman ρ={spearman.statistic:.3f}). Both top-10 lists include *Bradyrhizobium* and *Nocardioides*.

Nevertheless, **dominance and rank agreement are limited**: only {intersect[10]} genera occur in both top-10 sets, and only {intersect[20]} occur in both top-20 sets. Published CK-1 is led by several provisional labels (such as `PSRF01`, `SCGC-AG-212-J23` and `Gp6-AA40`) that have no exact-label counterpart in the PlusPF-8 Bracken output. Conversely, our reconstructed profile ranks *Streptomyces* and *Escherichia* much higher, and includes labels such as *Terribium* and *Mycolicibacterium* that are absent under those exact names in the published matrix. Distinct taxonomic naming frameworks and reference coverage likely contribute to this mismatch; an authoritative synonym/TaxID crosswalk would be needed to disentangle nomenclature from true assignment differences.

The independent Kraken2 step classified **1,833,747 / 33,808,912 non-host pairs (5.42%)** using the capped PlusPF-8 database. The known coverage limitation is described once here, rather than repeated throughout the biological interpretation. The attempt to use a larger, verified PlusPF-16 database was manually stopped for compute feasibility without producing a result; it cannot establish what PlusPF-16 would have classified. The authors' exact Kraken database build/size and percentage denominator are not stated, so disagreement is not evidence of an erroneous reconstruction, and the two genus profiles are not necessarily directly calibrated.

### One small Bracken accounting detail

The output contains `sum(new_est_reads) = {Btotal:,}`. The session's Bracken summary reports `1,830,146` genus-level reads before the final exported integer counts are aggregated, a difference of {1830146-Btotal:,} fragments (~{100*(1830146-Btotal)/1830146:.3f}%). Inspection of the installed Bracken abundance-estimation code confirms that redistributed abundance is calculated as floating point and each row is converted to an integer when written. Comparisons here renormalize those exported integer counts; `fraction_total_reads` is preserved verbatim alongside them for auditing.

## Figures

![Published top 20](fig1_published_top20.png)

![Bracken top 20](fig2_bracken_top20.png)

![Matched genera](fig3_shared_genera_comparison.png)

![Matched positive abundance scatter](fig4_shared_positive_scatter.png)

![Top-K overlap](fig5_topk_overlap.png)

## Reuse and reporting boundaries

- These are **single-sample**, profile-to-profile descriptive results. They do not support CK versus LS/MS significance tests.
- Unmatched names should be called **label differences**, not taxa proven absent in the other sample/profile.
- The abundance denominator and incomplete Kraken classification must appear in methods and limitations, but the Results section can straightforwardly discuss recovered taxa and shared ecological signals.
- Follow-up: inspect taxonomic synonyms and classification evidence for striking mismatches; potentially compare against a better-resourced reference on HPC.
- For the eventual Snakemake workflow, use the validated sample manifest, database paths/versions, independent source-matrix branch, and exact-label comparison rule shown here.

## Generated files

- `comparison_metrics.json`: summary statistics and explicit denominators.
- `all_genus_exact_label_comparison.csv`: all {len(names):,} union labels, ranks, normalized proportions and matching status.
- `shared_positive_genera.csv`: positive-abundance exact-name matches.
- `top20_each_profile.csv`: both top-20 lists with counterparts if present.
- `compare_ck1_profiles.py`: reproducible comparison from the extracted published CK-1 column and the Bracken output.
'''
    (out/'CK1_comparison_report.md').write_text(txt,encoding='utf-8')
    print('Saved to',out,'report/plots/comparison files. Shared positive:',len(both_pos),'spearman',spearman.statistic)

if __name__=='__main__':
    main()
