# differential-expression
Getting DESeq2 right: Experimental Design and the Stats that Power Differential Expression.

In this session, we will walk you through intermediate RNA-seq differential expression analysis in R using the popular DESeq2 package. You will learn the statistical foundations of experimental design and model specification in DESeq2 so you can obtain and interpret results that are robust and biologically meaningful for your study. Participants should be comfortable with R basics (e.g., working with data frames and factors). We will be using Positron to run our code.

You will learn how to:

- Connect research question with metadata file
- Assemble a DESeq2 dataset
- Design and explore a model matrix
- Perform basic data visualization (e.g. PCA plot, MA plot)
- Save and organize processed data files and results

## Requirements

- Basic knowledge of R programming language
- You may use your favorite IDE but if you would like to use Positron please be sure to have it installed before the class ([download from here](https://positron.posit.co/download.html)).

## Getting started

1. Download this repository: **Code → Download ZIP**, then unzip it.
2. In Positron, go to **File → Open Folder...** and select the `differential-expression` folder.
3. Install the packages (run once in the Console):

   ```r
   install.packages(c("BiocManager", "here"))
   BiocManager::install("DESeq2")
   ```

   DESeq2 is a Bioconductor package, so `install.packages("DESeq2")` will not find it.
4. Open `scripts/DESeq2.R` and run it from the top.

## Contents

```
differential-expression/
├── data/
│   ├── raw/                              # Data as downloaded from GREIN (GSE113754)
│   │   ├── metadata_Shank3.csv           # Sample metadata: genotype, condition
│   │   └── rawCounts_Shank3.csv          # Gene-level raw counts (genes × samples)
│   └── clean/                            # Cleaned versions of the raw data
│       ├── Shank3_metadata_clean.csv
│       └── Shank3_rawCounts_clean.csv
├── scripts/
│   └── DESeq2.R                          # Workshop script: DESeq2 analysis of the Shank3 data
├── presentation/
│   ├── differential-expression.qmd       # Slide source (Quarto, Beamer)
│   ├── differential-expression.pdf       # Rendered slides
│   ├── beamer-header.tex                 # Slide theme
│   └── assets/                           # Slide images
├── LICENSE
└── README.md
```

Running `scripts/DESeq2.R` creates a `results/` folder with the DESeq2 results tables and normalized counts.
