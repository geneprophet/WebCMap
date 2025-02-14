<img src="man/figures/overview.png" width="100%" />

# WebCMap

WebCMap is a highly efficient platform designed for query extended CMap,
renowned for its ultra-fast computation speed and minimal resource
consumption. It aims to provide users with a convenient and reliable
solution for connectivity mapping analysis. Its main features and
objectives include:

1.  Blazing-Fast Computation: Leveraging advanced algorithm optimization
    and an efficient web-accelerated CMap architecture, WebCMap achieves
    fast in connectivity mapping large-scale signature data,
    significantly enhancing research efficiency and saving users’ time.

2.  Low Resource Requirements: The platform is designed with a focus on
    optimizing computational resources, ensuring smooth operation even
    on standard hardware configurations, thereby lowering user barriers
    and hardware costs.

3.  Accurate Connectivity Analysis: Users can input query signature and
    efficiently mapping them with extended CMap. Six kinds of
    connectivity methods and another meta score are computed to
    accurately identify candidate signatures that are positively or
    negatively correlated with the query signature.

4.  Diverse Application Scenarios: The platform is suitable for various
    fields, including drug repurposing, disease mechanism research, and
    biomarker discovery, helping researchers quickly screen potential
    drug candidates and similar drug/toxic/compound searching.

## Installation

You can install the WebCMap like so:

``` R
devtools::install_github("geneprophet/WebCMap")
```

or download the source package

``` R
install.packages("/path/to/WebCMap_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Example

This is a basic example demonstrating how to screen candidate drugs for
coronary artery disease (CAD) using WebCMap:

``` R
# load WebCMap
library(WebCMap)
# load the example query_signature, which is derived from TWAS analysis of CAD (PMID:28209224)
data(query_signature)

head(query_signature)
#>   gene_id   zscore
#> 1    1723 3.321173
#> 2   51255 3.289684
#> 3   79609 3.153365
#> 4  122704 3.054415
#> 5  256126 2.839742
#> 6  114971 2.744463

# run the negative connectivity analysis
res = run_negative_CMap(query_signature = query_signature,K=50,cores=10)

head(res)
#>    signature_index Meta_score       WTCS      ES_up   ES_down  ES_up_padj
#> 17          198808          3 -0.4651067 -0.4791047 0.4511087 0.046995193
#> 18          211887          2 -0.4988176 -0.5193789 0.4782563 0.003771028
#> 34          519170          2 -0.5730391 -0.5524028 0.5936754 0.007146438
#> 2           629210          1 -0.5661617 -0.4765579 0.6557655 0.034321631
#> 4            32995          1 -0.4700981 -0.3736871 0.5665091 0.269430052
#> 16          625865          1 -0.5545614 -0.6064081 0.5027146 0.001681524
#>    ES_down_padj       XSum XSum_pvalue        CSS CSS_pvalue
#> 17 0.0481927711  -6.688600      0.0384 -0.1265102     0.0303
#> 18 0.0286858664  -8.044850      0.0248 -0.1428839     0.0133
#> 34 0.0071464378  -9.372175      0.0092 -0.1463990     0.0095
#> 2  0.0001320236 -12.208650      0.0026 -0.1340074     0.0215
#> 4  0.0033867776 -18.204400      0.0004 -0.1248601     0.0286
#> 16 0.0206405936  -7.201350      0.0099 -0.1609846     0.0053
#>    Spearman_correlation Pearson_correlation Consine_similarity       pert_id
#> 17           -0.2311988          -0.2524667         -0.2540377 BRD-K58660339
#> 18           -0.2027113          -0.2822857         -0.2778821 BRD-K86112430
#> 34           -0.1477813          -0.1996027         -0.1998024 BRD-A05334458
#> 2            -0.1703520          -0.2103692         -0.2118160 BRD-K31342827
#> 4            -0.1190621          -0.1426185         -0.1413339 BRD-K63750851
#> 16           -0.1761044          -0.1937879         -0.1923386 BRD-K92731339
#........
```

Visualization of CAD-lisofylline pair:

``` r
#filter the CAD-lisofylline pair:
data = res[res$cmap_name=="lisofylline",]
```

Radar plot:  

```R
radar_plot(data)
```





<img src="man/figures/README-Radar-1.png" width="100%" />

GSEA plot: 

```R
gsea_plot(data=data,query_signature=query_signature, K = 50)
```

<img src="man/figures/README-GSEA-1.png" width="100%" /> Venn plot: 

```R
venn_plot(data=data,query_signature=query_signature, K = 50)
```

 <img src="man/figures/README-Venn-1.png" width="100%" />
