#!/usr/bin/env python3
"""Reproduce the CK-1 exact-label comparison with Lv et al. (2026).

Requirements: Python 3, numpy, scipy.
Usage from the project root: python analysis/raw_read_comparison/compare_ck1_profiles.py
    --published data/processed/genus_abundance.tsv
    --bracken results/taxonomy/raw-read-pilot/CK1_pluspf8_c0_genus.bracken
    --outdir output
This script does not synonym-normalize taxa or infer taxon equivalence.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import pearsonr, spearmanr
from scipy.spatial.distance import braycurtis, jensenshannon


def write_csv(path, fields, rows):
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--published', type=Path, required=True)
    ap.add_argument('--bracken', type=Path, required=True)
    ap.add_argument('--outdir', type=Path, required=True)
    args = ap.parse_args()
    out = args.outdir
    out.mkdir(exist_ok=True, parents=True)
    with args.published.open(newline='', encoding='utf-8') as f:
        delimiter = '\t' if args.published.suffix == '.tsv' else ','
        p_rows = list(csv.DictReader(f, delimiter=delimiter))
    if delimiter == '\t':
        published = {r['genus'].strip(): float(r['CK-1']) for r in p_rows}
    else:
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

    print('Saved comparison tables and metrics to', out)


if __name__ == '__main__':
    main()
