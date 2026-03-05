# Transcriptome-scale Cas13 Knock-down (TraCK) libraries

## Environment

```
conda env create -f environment.yml
conda activate track
```

## Download expression datasets from DepMap and annotations from GENCODE

[DepMap Public 25Q3](https://depmap.org/portal/data_page/?tab=currentRelease)
```{bash}
mkdir -p datasets/expression

Download the following files into datasets/expression/:
- https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2025Q3&filename=OmicsExpressionTranscriptTPMLogp1HumanAllGenesStranded.csv
- https://depmap.org/portal/data_page/?tab=allData&releasename=DepMap%20Public%2025Q3&filename=Model.csv

python/convert_to_parquet.py --csv datasets/expression/OmicsExpressionTranscriptTPMLogp1HumanAllGenesStranded_25Q3.csv
```

## Run pipeline

First, configure according to your setup in `nextflow.config`.


## Visualization

module load nodejs/20.13.1-GCCcore-13.3.0
cd app
npm run dev

npm run preprocess
