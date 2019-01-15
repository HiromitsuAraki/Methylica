# Methylica
## GUI tool for independent component analysis of methylome data
Methylica is a GUI-based tool for independent component analysis (ICA) of methylome data generated by either whole-genome bisulfite sequencing (WGBS) or Infinium methylation array. Methylica not only provides ICA-based sample clustering but identifies independent components (ICs), or methylomic signatures, specific to various subsets of the samples. Major contributors to a subset-specific IC serve as methylation markers of the subset and imply biological processes underlying the IC. Methylica would thus be a powerful tool to analyze samples composed of multiple subtypes, such as those of cancers.
<br>

## Install/Launch Methylica
1.  Install [R environment](https://www.r-project.org/)
2.  Install [shiny](https://shiny.rstudio.com).  
`install.packages("shiny")`
3.  Launch Methylica  
The following R code will launch Methylica.  
`shiny::runGitHub("HiromitsuAraki/Methylica")`
<br>

## Input file format
- Methylome data
  - 1st column: Chr
  - 2nd column: Start
  - 3rd column: End
  - 4th column: Gene symbol
  - 5th column ~ : Methylome data of each sample
  <br>
- Sample meta data
  - 1st column: Sample ID
  - 2nd column ~ : Status of the features (e.g. cancer subtype, stage, gender)  
  **NOTE: The status of the features should be discrete, as Methylica cannot accept metadata with continuous values (e.g. age, tumor size, and survival date).**  
<br>

## Implementations
### Data uploading
Methylica requires methylome data and sample metadata as its inputs. The former is a matrix of methylation levels, rows and columns of which correspond to genomic regions and samples, respectively. The latter is a tab-delimited text file, rows and columns of which correspond to samples and features (e.g. gender, risk factor, and cancer subtype), respectively. The status of the features should be discrete, as Methylica cannot accept metadata with continuous values (e.g. age, tumor size, and survival date).

<img src="./README_files/Figures/DataUpload.png" width=300x300>
<br>

### Parameter setting
Methylica requires methylome data and sample metadata as its inputs. The former is a matrix of methylation levels, rows and columns of which correspond to genomic regions and samples, respectively. The latter is a tab-delimited text file, rows and columns of which correspond to samples and features (e.g. gender, risk factor, and cancer subtype), respectively. The status of the features should be discrete, as Methylica cannot accept metadata with continuous values (e.g. age, tumor size, and survival date).

<img src="./README_files/Figures/Parameters.png" width=300x300>
<br>

### Visualization ICs
#### Heatmap clustering of ICs
Methylica requires methylome data and sample metadata as its inputs. The former is a matrix of methylation levels, rows and columns of which correspond to genomic regions and samples, respectively. The latter is a tab-delimited text file, rows and columns of which correspond to samples and features (e.g. gender, risk factor, and cancer subtype), respectively. The status of the features should be discrete, as Methylica cannot accept metadata with continuous values (e.g. age, tumor size, and survival date).

<img src="./README_files/Figures/HeatmapClustering.png" width=500x500>
<br>

#### Boxplots
Methylica requires methylome data and sample metadata as its inputs. The former is a matrix of methylation levels, rows and columns of which correspond to genomic regions and samples, respectively. The latter is a tab-delimited text file, rows and columns of which correspond to samples and features (e.g. gender, risk factor, and cancer subtype), respectively. The status of the features should be discrete, as Methylica cannot accept metadata with continuous values (e.g. age, tumor size, and survival date).

<img src="./README_files/Figures/Boxplots.png" width=500x500>
<br>

#### Loading factor table
Methylica requires methylome data and sample metadata as its inputs. The former is a matrix of methylation levels, rows and columns of which correspond to genomic regions and samples, respectively. The latter is a tab-delimited text file, rows and columns of which correspond to samples and features (e.g. gender, risk factor, and cancer subtype), respectively. The status of the features should be discrete, as Methylica cannot accept metadata with continuous values (e.g. age, tumor size, and survival date).

<img src="./README_files/Figures/ICtable.png" width=500x500>
<br>


### Highliting highly contributed regions
#### Heatmap clustering of highly contributed regions
Methylica requires methylome data and sample metadata as its inputs. The former is a matrix of methylation levels, rows and columns of which correspond to genomic regions and samples, respectively. The latter is a tab-delimited text file, rows and columns of which correspond to samples and features (e.g. gender, risk factor, and cancer subtype), respectively. The status of the features should be discrete, as Methylica cannot accept metadata with continuous values (e.g. age, tumor size, and survival date).

<img src="./README_files/Figures/HeatmapClustering_highLF.png" width=500x500>
<br>

#### Loading factor table
Methylica requires methylome data and sample metadata as its inputs. The former is a matrix of methylation levels, rows and columns of which correspond to genomic regions and samples, respectively. The latter is a tab-delimited text file, rows and columns of which correspond to samples and features (e.g. gender, risk factor, and cancer subtype), respectively. The status of the features should be discrete, as Methylica cannot accept metadata with continuous values (e.g. age, tumor size, and survival date).

<img src="./README_files/Figures/ICtable_highLF" width=500x500>
<br>

#### Boxplots
Methylica requires methylome data and sample metadata as its inputs. The former is a matrix of methylation levels, rows and columns of which correspond to genomic regions and samples, respectively. The latter is a tab-delimited text file, rows and columns of which correspond to samples and features (e.g. gender, risk factor, and cancer subtype), respectively. The status of the features should be discrete, as Methylica cannot accept metadata with continuous values (e.g. age, tumor size, and survival date).

<img src="./README_files/Figures/Boxplots.png" width=500x500>
<br>
