seq_autoqc
==========

These routines run on data produced by bamcheck and are used as part of the autoqc process in HGI pipelines.  

bamcheck files used to be produced by the `bamcheck` program but are now incorporated into `samtools` and produced by the `samtools stats` command (as of the 1.0 release candidates beginning with "0.2.0-rc1").  

bamcheckr
---------
An R package for reading and writing bamcheck data, and for calculating indel_peaks, 
quality_dropoff, and base_content_deviation.

This package can be installed directly from git using `devtools`
```R
library(devtools)
install_github("samstudio8/seq_autoqc", subdir="bamcheckr")
```

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

If you encounter trouble (`could not find function "loadMethod"`) attempting use with RScript you could try the somewhat archaic `R CMD BATCH` as follows:
```bash
R CMD BATCH '--args input.bamcheck output.bamcheck --plot-base-path=plots/bamplot' bamcheck_augment_summary.R
```
