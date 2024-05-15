
# SigRM

<!-- badges: start -->
<!-- badges: end -->

The SigRM package, a novel statistical framework for the effective mining of single-cell epitranscriptomics datasets comprised of a very large number of cells.   Specifically designed for discerning m6A modification sites, SigRM relies on four independent negative binomial distributions, with linked variances and means through local regressions.   This comprehensive toolkit comprises three pivotal functions.   Firstly, SigRMtest serves as the cornerstone for analyzing single-cell RNA methylation sites, providing a deep dive into the intricacies of m6A modifications within individual cells.   Secondly, SigRM_similarity_test assumes a central role in this analysis, enabling the selection of control groups based on the similarity of each test cell.   Users have the flexibility to customize the number of control cells.   Lastly, SigRM_cluster_test is essential for clustering cells based on gene expression or read counts, accommodating both test and control groups seamlessly.   Effective utilization of this function necessitates prior data clustering, with each cluster ideally embodying a mix of test and control cells, enabling a comprehensive assessment of test cell behavior.   This innovative approach relies on control cells with similar gene expression patterns as a baseline for comparison, systematically evaluating test cell behavior within their respective clusters.   SigRM not only achieved improved performance in RNA methylation site detection compared to generic and state-of-the-art models, but also provided a rigorous quantification of the RNA methylation level that facilitate various forms of downstream analysis.  For any inquiries, please do not hesitate to contact Haozhe.Wang17@student.xjtlu.edu.cn.


## Installation

You can install the development version of SigRM from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("whz991026/SigRM")
```

## data

The files associated with this study include:

*SNP File (SNP_all.rds):

Data Details: Stored in RDS format, containing a GRanges object.
Content: Each entry in the GRanges object represents a specific genomic region corresponding to an m6A modification site detected through scDART-seq.
seqnames: Factor Rle object containing chromosome or genomic sequence names.
ranges: IRanges object containing genomic intervals (start and end positions).
strand: Factor Rle object indicating the strand (directionality) of the genomic region.
mcols: DataFrame object containing optional metadata columns, such as quality scores, coverage depth, and mutation data.
seqinfo: Seqinfo object providing information about the genomic sequences present in the GRanges object.
Purpose: Provides detailed genomic information about identified m6A modification sites, facilitating further analysis of their distribution, characteristics, and genomic context in HEK293T cells.


*Frequency File (frequency_all.rds):

Data Details: Also stored in RDS format, consisting of a list.
Content: Each item in the list represents a single-cell, containing counts of methylated and unmethylated reads for corresponding m6A modification sites detected in scDART-seq data.
Purpose: Offers quantitative data on the abundance of methylated and unmethylated sequences at each m6A site across individual cells, enabling investigation of m6A modification patterns at a single-cell level.


*Expression TPM File (expression_TPM.rds):

Data Details: Stored as an RDS file, comprising a list.
Content: Each item in the list represents a single-cell, with corresponding TPM values for gene expression.
Purpose: Provides information on gene expression levels across individual cells, facilitating examination of potential correlations between m6A modification patterns and gene expression profiles in scDART-seq data from HEK293T cells.

*Gene Information File (gene_informations.rds):

Data Details: Stored as an RDS file, comprising a data frame.
Content: Includes information such as Gene ID, Gene Name, Reference, Strand, Start position, End position, and Coverage.
Purpose: Offers additional details about gene expression data, aiding in the interpretation and analysis of gene expression profiles in conjunction with m6A modification patterns.




All three files can access on: https://doi.org/10.6084/m9.figshare.25815862

with DOI: 10.6084/m9.figshare.25815862



## Example

It demonstrates how to utilize the main function, SigRMtest(), from the SigRM package. The primary input for the SigRMtest function is a matrix or data frame containing four sets of read counts, where rows represent genes, and columns represent single cells. These four sets of read counts are ordered as follows: methylation read counts in the control group, methylation read counts in the test group, unmethylation read counts in the control group, and unmethylation read counts in the test group. This data can be obtained by calling the simulateData() function.

The output of the SigRMtest function is a list with a length of 7:

* The first part consists of a data frame representing methylation proportions (with rows for genes and columns for single cells) in the test cells.
* The second part comprises a mean methylation proportion vector (with rows for genes and columns for single cells) in the control cells.
* The third part contains a data frame displaying log2 risk ratios for the test cells.
* The fourth part includes a data frame illustrating log2 odds ratios for the test cells.
* The fifth part encompasses a data frame containing p-values for the test cells.
* The sixth part represents estimated gene abundances.
* The seventh part is a data frame of adjusted p-values for the test cells. 


```{r}
set.seed(1)
# first simulate the data
data <- simulateData(test_num = 10,control_num = 30)
# put into the main function SigRMtest
res <- SigRMtest(data[[1]],data[[2]],data[[3]],data[[4]])

# meth proportion data frame (row is gene and column is single cell) in the test cells. 
head(res[[1]])

# mean meth proportion vector (row is gene and column is single cell) in the control cells.
head(res[[2]])

# log2 risk ratio data frame for the test cells. 
head(res[[3]])

# log2 odds ratio data frame for the test cells. 
head(res[[4]])

# p value data frame for the test cells. 
head(res[[5]])

# estimated gene abundance
head(res[[6]])

# adjusted p value data frame for the test cells
head(res[[7]])
```
