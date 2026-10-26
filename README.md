# Clinical AI & Computational Medicine

Runnable, **synthetic-data** computational medicine projects demonstrating biomedical data engineering, numerical modeling, and reproducible validation.

| Project | Implemented methods | Validation |
|---|---|---|
| [Clinical genomics](clinical-genomics/) | Pedigree inference, CNV grouping, cis/trans constraint graph, haplotype enumeration, genetic diagnosis classification | Full reference agreement on 4,944 alleles and 838 diagnosis rows |
| [Vancomycin pharmacokinetics](vancomycin-pk/) | Infusion reconstruction, level revision filtering, two-compartment ODE transitions, constrained PK parameter fitting, AUC prediction, regimen ranking | PK fitting/observed AUC validated; **projected regimens still differ from reference** |

Each project includes a readable analysis script (`analyze.py`), full computational pipeline (`full_pipeline.py`), synthetic example and full-size datasets, tests, README, and generated results. See project READMEs for instructions and limitations.

This portfolio does not represent a clinically certified system. **Do not use it for diagnosis, prescribing, or patient care.**

## Additional past work

Comparative anticoagulant bleeding risk and sepsis treatment-timing analysis are described elsewhere and not represented by executable code in this repository; no implementation is claimed here.

## Publication rights

The synthetic examples were prepared in connection with technical evaluation work. Confidential evaluator outputs are excluded.
