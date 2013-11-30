seq_autoqc
==========

These routines run on data produced by bamcheck and are used as part of the autoqc process in HGI pipelines.


bamcheckr
---------
An R package for reading and writing bamcheck data, and for calculating indel_peaks, 
quality_dropoff, and base_content_deviation.

bamcheck_augment_summary.R 
--------------------------
An R script that takes bamcheck as input and produces augmented bamcheck as output, in which 
additional Summary Number (SN) entries are added for indel peaks, quality dropoff, and base 
content deviation. 

You can obtain help by running: 
```bash
Rscript bamcheck_augment_summary.R --help
```
 
Usage with defaults is as simple as:
```bash
Rscript bamcheck_augment_summary.R input.bamcheck output.bamcheck
```

These defaults are equivalent to:
```bash
Rscript bamcheck_augment_summary.R --indel-runmed-k=25 --indel-baseline-method=runmed --base-content-runmed-k=25 --base-content-baseline-method=mean --quality-dropoff-runmed-k=25 --quality-dropoff-ignore-edge-cycles=3 --quality-dropoff-high-iqr-threshold=1 input.bamcheck output.bamcheck
```
