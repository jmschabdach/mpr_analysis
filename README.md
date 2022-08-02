# About

Jenna Young (@jmschabdach)

Take FreeSurfer outputs and convert them into a format that can be used more easily in Python/R analyses.

# Requirements

In Python

```
pandas
```

In R

- `ggplot2`
- `patchwork`

# Use

## Step 1: Preprocess the MPRAGEs

Make sure the MPRAGE images are in an optimal format for use with FreeSurfer's reconall tool.

For a single image: 

`bash preprocWashUACPCAlignment.sh [args]`

For a set of images on a cluster:

` `

## Step 2: Run the MPRAGEs through FreeSurfer

Run `recon-all` 

## Step 3: Grab the statistics produced by FreeSurfer and put them in a database file

Run the following script to extract stats measures from the FreeSurfer output of each image in the path:

`extractFreeSurferMeasures.py -p /path/to/freesurfer/output/directory`

## Step 4: Use the .db file to produce graphs describing the data

Use `queryFreeSurferMeasures.py`

# Statistical Analysis

1. Combine the regular and infant FS outputs using `joinFsIfsDataframes.py`
2. In R, modify the filenames etc and run `load_clip_dataframe.r` to filter and format the dataframe
3. In R, modify the filenames etc and run `run_combat.r` and `runNeuroHarmonize.py` as comments indicate to combat the data
4. In R, modify the filenames and title and run `create_scatterplots.r` to generate basic scatterplots of selected phenotypes colored by primary reason for scan.
