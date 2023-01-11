
# DWLS

<!-- badges: start -->
<!-- badges: end -->

Forked from https://bitbucket.org/yuanlab/dwls/src/master/ and modified for speed.

## Installation

You can install the development version of DWLS like so:

``` r
remotes::install_github('alexvpickering/dwls')
```

## Example

DWLS requires 3 inputs:
- cluster labels and marker genes
- single-cell counts
- bulk counts

Marker genes in the appropriate format can be obtained from single-cell counts in a variety of ways. For example

``` r
library(DWLS)

library(Seurat)
data(pbmc_small)
id <- Idents(pbmc_small)
scdata <- pbmc_small[['RNA']]@counts

markers <- DEAnalysisPresto(scdata, id)
```

The next step is to build the single-cell signature matrix:

```r
sig <- buildSignatureMatrix(scdata, id, markers)
```

Next run deconvolution:

```r
# example bulk signature
bulk <- rowSums(scdata)

# trim to common genes
tr <- trimData(sig, bulk)

# run deconvolution
res <- solveDampenedWLS(tr$sig, tr$bulk)

# parallel implementation for multiple bulk signatures
bulk2 <- cbind(bulk, bulk)
tr2 <- trimData(sig, bulk2)
res2 <- solveParallel(tr2$sig, tr2$bulk)
```

