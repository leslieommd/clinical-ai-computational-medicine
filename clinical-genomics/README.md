# Clinical Genomics: Pedigree-Aware Haplotype Reconstruction

A fully runnable Python implementation of autosomal-recessive allele reconstruction from synthetic pedigree, sequence-variant, deletion, and phasing evidence.

## Methods

1. Group CNV records using dosage and reciprocal genomic overlap; preserve raw identifiers.
2. Constrain inheritance with parental observations and homozygous events.
3. Model cis/trans constraints as a parity graph and detect contradictions.
4. Enumerate allowable orientations and classify each patient-gene combination across all compatible assignments.
5. Export allele-level assignments, diagnoses, and QC counts.

## Run

```bash
python full_pipeline.py --input data_full --output results_full
python -m unittest discover -s tests -v
```

The smaller `data/` directory provides a quick example. Full synthetic data are under `data_full/`.

**Full synthetic run:** 4,944 event assignments; 838 patient-gene diagnoses; 72 phase-conflict groups; 118 merged CNV observations; 104 unresolved diagnoses. Results have been checked against a separately held reference. The data do not correspond to clinical patients.

This is a simplified synthetic Mendelian model, not clinical variant interpretation following ACMG/AMP guidelines, and it must not be used for patient diagnosis.
