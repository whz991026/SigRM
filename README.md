
# SIGMR

<!-- badges: start -->
<!-- badges: end -->

The SIGMR package presents a novel statistical approach, utilizing a 
    similarity-guided multi-resolution analysis for the detection of modification 
    sites in DART-seq data. Specifically designed for discerning m6A modification 
    sites, SIGMR relies on four independent negative binomial distributions, with 
    linked variances and means through local regressions. This comprehensive toolkit 
    comprises three pivotal functions. Firstly, SIGMRtest serves as the cornerstone 
    for analyzing single-cell RNA methylation sites, providing a deep dive into the 
    intricacies of m6A modifications within individual cells. Secondly, 
    SIGMR_similarity_test assumes a central role in this analysis, enabling the 
    selection of control groups based on the similarity of each test cell. Users 
    have the flexibility to customize the number of control cells. Lastly, 
    SIGMR_cluster_test is essential for clustering cells based on gene expression 
    or read counts, accommodating both test and control groups seamlessly. 
    Effective utilization of this function necessitates prior data clustering, 
    with each cluster ideally embodying a mix of test and control cells, enabling 
    a comprehensive assessment of test cell behavior. This innovative approach 
    relies on control cells with similar gene expression patterns as a baseline 
    for comparison, systematically evaluating test cell behavior within their 
    respective clusters. The SIGMR package equips researchers with powerful tools 
    to enhance the accuracy of methylation site results derived from DART-seq data. 
    For any inquiries, please do not hesitate to contact Haozhe.Wang17@student.xjtlu.edu.cn.

## Installation

You can install the development version of SIGMR from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("whz991026/SIGMR")
```

## Example

It demonstrates how to utilize the main function, SIGMRtest(), from the SIGMR package. The primary input for the SIGMRtest function is a matrix or data frame containing four sets of read counts, where rows represent genes, and columns represent single cells. These four sets of read counts are ordered as follows: methylation read counts in the control group, methylation read counts in the test group, unmethylation read counts in the control group, and unmethylation read counts in the test group. This data can be obtained by calling the simulateData() function.

The output of the SIGMRtest function is a list with a length of 7:

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
# put into the main function SIGMRtest
res <- SIGMRtest(data[[1]],data[[2]],data[[3]],data[[4]])

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
