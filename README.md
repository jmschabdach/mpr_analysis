# About

Jenna Young (@jmschabdach)

Develop brain growth charts using FreeSurfer and SynthSeg phenotype outputs.

# Requirements

In Python

```
pandas
```

In R

- `ggplot2`
- `patchwork`
- `gamlss`

# Use

There are two uses for this repo. The first is to completely replicate the results presented in the original SLIP paper. The second is to use the functions written for the original SLIP paper analysis to build your own growth charts. We present the second use case first as it is likely of more interest to the majority of folks visiting this page. 

# Build Your Own Growth Charts

Uses `./r/build_your_own_growth_chart.R` and `./r/lib_mpr_analysis.r` 


The `build_your_own_growth_chart.R` script was assembled to walk a user through building a GAMLSS model using their own neuroimaging phenotype and demographic data. If you customize it and use it in your research, we politely ask that you cite the following papers (at least the first one if your citations are limited):

[SLIP paper](https://www.medrxiv.org/content/10.1101/2023.01.13.23284533v1)

[Lifespan Nature paper](https://www.nature.com/articles/s41586-022-04554-y)


The models build in this script require the following data:
- logAge: numeric type with log(post conception age in days, base=10). In (1), we use a conversion factor of 325.25 days/year and a post conception offset of 280 days if post conception age is not available.
- sex: factor with M or F values.
- phenotype: numeric type with the measurements of the phenotype of interest

For FreeSurfer (FS) neuroimaging phenotypes, the measurement SurfaceHoles is included in the GAMLSS model. SynthSeg (SS) does not produce SurfaceHoles.

The script is commented to walk the user through the following steps:
1. Loading data from a .csv file
2. Building a growth chart GAMLSS model
3. Estimating the median centile for the specified phenotype
4. Estimating the age at the peak value of the median centile of the specified phenotype
5. Estimating a set of centile curves from the GAMLSS file
6. Plotting the original data and centile curves
7. Estimate the centile corresponding to each phenotype value in the data
8. Plotting the phenotype centiles on a violin plot

# SLIP Paper Analysis

This repo assumes that the anatomical MRIs to be included in the brain growth charts have been processed using FreeSurfer 6.0/Infant FreeSurfer (with an age cutoff of 2 years old) and SynthSeg+ independently. 

### Phenotype File Building

If the directories mirror the following format, the path variables in the scripts contained within the `extractPhenotypes` directory can be modified and used to extract all phenotypes from FreeSurfer/Infant FreeSurfer and SynthSeg output directories. 

```
data/
|-- derivatives/
    |-- freesurfer/
        |-- sub-001/
        |-- sub-002/
    |-- sythseg/
        |-- sub-001/
        |-- sub-002/

```

Run the following script to extract stats measures from the FreeSurfer output of each image in the path:

`extractFreeSurferMeasures.py -p /path/to/data/derivatives/freesurfer`

The phenotypes extracted in this step contain local and global phenotypes for each subject as well as demographic information about the subjects.

### Statistical Analysis

Run the following files in order:
```
01_load_clip_dataframe.r
02_run_combat.r
03_create_centile_correlation_plots.r
03b_run_demographic_analysis.r
03c_calculate_grader_agreement.r
```

The files beginning with `03` can be run in parallel: they all use the outputs from step 02 (which uses outputs from step 01) and are independent from each other.

The `04_create_outofsample_ss_plots.R` file was used to evaluate the fit of the SynthSeg models to out of sample data not included in the original GAMLSS models built in `03_create_centile_correlation_plots.r`.
