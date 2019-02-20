# Methylica
## GUI tool for independent component analysis of methylome data
*Methylica* is a GUI-based tool for independent component analysis (ICA) of methylome data generated by either sequence-based methylome data, such as whole-genome bisulfite sequencing (WGBS) or Reduced Representation Bisulfite Sequencing (RRBS), or Infinium methylation array. *Methylica* not only provides ICA-based sample clustering but identifies independent components (ICs), or methylomic signatures, specific to various subsets of the samples. Major contributors to a subset-specific IC serve as methylation markers of the subset and imply biological processes underlying the IC. *Methylica* would thus be a powerful tool to analyze samples composed of multiple subtypes, such as those of cancers.
<br>

## Install/Launch Methylica
1.  Install [R environment](https://www.r-project.org/)
2.  Install [shiny](https://shiny.rstudio.com).  
`install.packages("shiny")`
3.  Launch *Methylica*  
The following R code will launch *Methylica*.  
`shiny::runGitHub("HiromitsuAraki/Methylica")`
<br>

## Input file format
*Methylica* requires methylome data and sample metadata as its inputs. The former is a matrix of methylation levels, rows and columns of which correspond to genomic regions and samples, respectively. The latter is a tab-delimited text file, rows and columns of which correspond to samples and features, respectively. The status of the features should be discrete, as *Methylica* cannot accept sample metadata with continuous values (e.g. survival date). We provide [*MatrixMaker*](https://github.com/HiromitsuAraki/MatrixMaker), which converts user's methylome data to a data matrix compatible with *Methylica*.
- Methylome data
  - 1st column: Chr
  - 2nd column: Start
  - 3rd column: End
  - 4th column: Gene symbol (if users analyize CpG island or CpG island shores, "CGI" should be described in 4th column)
  - 5th column ~ : Methylome data of each sample
  <br>
- Sample meta data
  - 1st column: Sample ID
  - 2nd column ~ : Status of the features (e.g. cancer subtype, stage, gender)  
  
**NOTE: The sample IDs on 5th column ~ in methylome data and on 1st column in sample meta data have to be identical.**  
<br>

## Implementations
### Data uploading
*Methylica* requires methylome data and sample metadata as its inputs. Please refer **Input file format** about the file format of methylome data and sample metadata for *Methylica*. Users need to assign the file location from a browser as below.

<img src="./README_files/Figures/DataUpload.png" width=300x300>
<br>

### Parameter setting
Following data upload, *Methylica* requests its users to select a platform (Sequence-based or Infinium methylation array), species with its reference genome version, genomic elements to be analyzed (CpG island, CpG island shores, gene body, first intron and promoter), and k or the number of ICs (minimum = 2; maximum = the number of samples). When users select CpG island shores, which are 2kbp flanking regions on each side of CpG island, genomic coordinates of input data matrix have to be identical to CpG islands. We recommend that data matrix to be analyzed is generated by *MatrixMaker*. *Methylica* provides a default setting of k, defined as the first k components whose cumulative contribution ratio exceeds 80% in principal component analysis. When users select all parameters, users need to press "Run" button to start analysis.

<img src="./README_files/Figures/Parameters.png" width=300x300>
<br>

### Visualizing ICs
When ICA is completed, *Methylica* provides a heatmap of clustering based on the mixing coefficient matrix A, boxplots to compare weightings among the user-defined sample features, and a table of loading values comprising each IC.

#### Heatmap clustering of mixing coefficient matrix
The vertical axis corresponds to mixing coefficients (sample weight in each IC) and the horizontal axis corresponds to samples. Heatmap clustering can be downloaded by pressing *Download* button.

<img src="./README_files/Figures/HeatmapClustering.png" width=500x500>
<br>

#### Boxplots of mixing coefficient matrix
The x axis is sample weights in each IC and the y axis is a sample feature defined in sample meta data. Users can select another sample feature in *Sample property* box as a sample feature defined in sample meta data box. Boxplots can be downloaded by pressing *Download* button.

<img src="./README_files/Figures/Boxplots.png" width=500x500>
<br>

#### Loadings table
Regions with absolute loadings >=2 are appeared in this table. Users can retrieve user's interesting genes in *Search* text box.

<img src="./README_files/Figures/ICtable.png" width=500x500>
<br>


### Highlighting highly contributed regions
Once an IC of interest is identified, *Methylica* can display its high loading genomic regions as a methylation heatmap and a clickable list to select a region whose methylation levels in sample features are indicated as boxplots. The list can be exported for various analyses (e.g. GO enrichment and sequence motif search) to facilitate biological interpretation of the IC.

#### Heatmap clustering of highly contributed regions
Users can assign IC and threshold of loadings to be visualized in "IC" and "Loadings" text box respectively. The vertical axis corresponds to genomic regions with high loadings and the horizontal axis corresponds to samples. The genomic regions are sorted by loadings. Heatmap clustering can be downloaded by pressing *Download* button.

<img src="./README_files/Figures/HeatmapClustering_highLF.png" width=500x500>
<br>

#### Loadings table of highly contributed regions
Boxplot of methylome data in each region is appeared by clicking gene symbol (see **Boxplot of highly contributed region**). Whole list of loadings table can be downloaded by pressing *Download* button.

<img src="./README_files/Figures/ICtable_highLF.png" width=500x500>
<br>

#### Boxplot of highly contributed region
The x axis is methylation level and the y axis is sample features. Users can select another sample feature defined in sample meta data in *Sample property* box. 

<img src="./README_files/Figures/Boxplot_highLF.png" width=500x500>
<br>
