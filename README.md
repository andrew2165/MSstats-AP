# MSstatsAP personal fork

The GitHub repository is `MSstats-AP`; the installed R package is
`MSstatsAP`. R package names cannot contain hyphens. The AP packages install
alongside the original packages without replacing them. Function names remain
unchanged; use explicit namespace calls when both versions are attached.

Install the renamed forks in dependency order from their feature branches
(after these changes have been committed and pushed):

```r
remotes::install_github("andrew2165/MSstats-AP@bait-normalization-release-3.23")
remotes::install_github("andrew2165/MSstatsTMT-AP@backport/bait-normalization-3.23")
remotes::install_github("andrew2165/MSstatsShiny-AP@feature/label-free-bait-normalization")

library(MSstatsAP)
library(MSstatsTMTAP)
library(MSstatsShinyAP)
MSstatsShinyAP::launch_MSstatsShiny()

# Check the actual installed package locations and GitHub source commits.
for (pkg in c("MSstatsAP", "MSstatsTMTAP", "MSstatsShinyAP")) {
  print(find.package(pkg))
  print(packageDescription(pkg)$RemoteSha)
}
```

Other MSstats family dependencies retain their original package names and may
also load the upstream MSstats namespace. AP analyses use the AP imports and
namespace calls. The upstream overview below is retained for background.
The linked upstream websites describe the original packages.

---

MSstats
=======

<!-- badges: start -->
  [![Travis build status](https://travis-ci.org/Vitek-Lab/MSstats-dev.svg?branch=master)](https://travis-ci.org/Vitek-Lab/MSstats-dev)
[![Codecov test coverage](https://codecov.io/gh/Vitek-Lab/MSstats-dev/branch/master/graph/badge.svg)](https://codecov.io/gh/Vitek-Lab/MSstats-dev?branch=master)
<!-- badges: end -->

MSstats is an R-based/Bioconductor package for statistical relative quantification of peptides and proteins in mass spectrometry-based proteomic experiments. It is applicable to multiple types of sample preparation, including label-free workflows, workflows that use stable isotope labeled reference proteins and peptides, and work-flows that use fractionation. It is applicable to targeted Selected Reactin Monitoring(SRM), Data-Dependent Acquisiton(DDA or shotgun), and Data-Independent Acquisition(DIA or SWATH-MS). This github page is for sharing source and testing. The official webpage is msstats.org
