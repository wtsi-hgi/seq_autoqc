seq_autoqc
==========

These routines run on data produced by bamcheck and are used as part of the autoqc process in HGI pipelines.


indel_vs_readcycle_peaks
------------------------
Issues that affect a particular read-cycle (such as temporary bubbles) tend to cause large spikes in the indel vs readcycle plots. This script detects those spikes.


quality_dropoff
---------------
Sometimes the quality of a run drops off early. This can be due to increasing loss of phasing as bases are incorporated or skipped when conditions are not ideal (e.g. because of bad reagents, fluidics issues, temperature, etc).

```bash
R --vanilla --slave --args bamcheck="${bamcheck_filename}" iqr.threshold=10 fail.thresh.ccc=12 warn.thresh.ccc=5 outdir="${outdir}" < pass.fail.warn.iqr.test.R
```
 
