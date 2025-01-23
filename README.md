
<!-- README.md is generated from README.Rmd. Please edit that file -->

# WebCMap

<!-- badges: start -->
<!-- badges: end -->

WebCMap is a highly efficient platform designed for query extended CMap,
renowned for its ultra-fast computation speed and minimal resource
consumption. It aims to provide users with a convenient and reliable
solution for connectivity map analysis. Its main features and objectives
include:

1.Blazing-Fast Computation: Leveraging advanced algorithm optimization
and an efficient computational architecture, WebCMap achieves fast
responses in comparing and analyzing large-scale signature data,
significantly enhancing research efficiency and saving users’ time.

2.Low Resource Requirements: The platform is designed with a focus on
optimizing computational resources, ensuring smooth operation even on
standard hardware configurations, thereby lowering user barriers and
hardware costs.

3.Accurate Connectivity Analysis: Users can upload query signatures
(query signatures) and efficiently compare them with extended CMap. Six
kinds of connectivity scores are computed to accurately identify
candidate signatures that are positively or negatively correlated with
the query signature.

4.Diverse Application Scenarios: The platform is suitable for various
fields, including drug repurposing, disease mechanism research, and
biomarker discovery, helping researchers quickly screen potential drug
candidates and similar drug/toxic/compound searching.

## Installation

You can install the development version of WebCMap like so:

``` r
# FILL THIS IN! HOW CAN PEOPLE INSTALL YOUR DEV PACKAGE?
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(WebCMap)
#> 
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
