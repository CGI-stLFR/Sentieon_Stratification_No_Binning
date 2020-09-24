# Stratification Pipeline

This pipeline was developed to handle a Sentieon assembly that wasn't phased.
That said they actually meant haploid and this pipeline was only used once and is kinda worthless.
I'm adding some comments just in case

## Running the pipeline

Symlink the `Assembly/` directory under `Assembly`.
This pipeline is primarily for simulated libraries.
Be sure to name it appropriately.
Add the symlinked directory name to the config file under `samples`.
Then just run snakemake as below.

```
# -s supplies the snakefile
# --configfile supplies the config file
# -j supplied threads
snakemake -s run_binning.Snakefile --configfile run_binning.config -j 60 2>&1 | tee snakemake.err.txt
```

## run_binning.config

- samples
    - samples specifies which samples to evaluate
    - just add the name of the dir that's been added to `Assembly` before running the pipeline
