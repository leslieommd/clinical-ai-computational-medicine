"""Exploratory quality-control summary for synthetic family genomics data.
Run: python analyze.py
"""
from collections import Counter, defaultdict
from pathlib import Path
import csv, json

DATA = Path(__file__).parent / 'data'

def load(name):
    with (DATA / name).open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f))

pedigree, variants, cnvs, phases, genes = [load(n) for n in ('pedigree.csv','variants.csv','cnvs.csv','phase.csv','genes.csv')]
people = {r['sample_id'] for r in pedigree}
probands = {r['sample_id'] for r in pedigree if r['mother_id'] and r['father_id']}
assert all(r['mother_id'] in people and r['father_id'] in people for r in pedigree if r['sample_id'] in probands)
assert all(r['sample_id'] in people for r in variants + cnvs + phases)
summary = {
    'families': len({r['family_id'] for r in pedigree}),
    'proband_count': len(probands),
    'variant_records': len(variants),
    'cnv_records': len(cnvs),
    'phase_records': len(phases),
    'genes': sorted({r['gene'] for r in genes}),
    'variant_classes': dict(sorted(Counter(r['classification'] for r in variants).items())),
    'phase_relationships': dict(sorted(Counter(r['relationship'] for r in phases).items())),
    'proband_gene_event_counts': len({(r['sample_id'],r['gene']) for r in variants+cnvs if r['sample_id'] in probands}),
}
print(json.dumps(summary,indent=2))
