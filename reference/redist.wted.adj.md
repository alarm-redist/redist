# Create Weighted Adjacency Data

Create Weighted Adjacency Data

## Usage

``` r
redist.wted.adj(map = NULL, plans = NULL)
```

## Arguments

- map:

  redist_map

- plans:

  redist_plans

## Value

tibble

## Examples

``` r
data(iowa)
shp <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.01)
plans <- redist_smc(shp, 100)
#> SEQUENTIAL MONTE CARLO
#> Sampling 100 99-unit maps with 4 districts and population between 753973 and 769205.
#> Split [0/3] ■                                | ETA?
#> Split [3/3] ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■  | ETA 0s
#> 
redist.wted.adj(shp, plans = plans)
#> Simple feature collection with 222 features and 3 fields
#> Geometry type: LINESTRING
#> Dimension:     XY
#> Bounding box:  xmin: 4192376 ymin: 2974708 xmax: 5729102 ymax: 3978011
#> Projected CRS: NAD83(HARN) / Iowa North (ftUS)
#> # A tibble: 222 × 4
#> # Rowwise: 
#>        i     j                           geometry    wt
#>  * <int> <int>      <LINESTRING [US_survey_foot]> <dbl>
#>  1     1     2 (4654501 3220689, 4590239 3111565)  0.94
#>  2     1    15 (4654501 3220689, 4529015 3222722)  0.96
#>  3     1    39 (4654501 3220689, 4647780 3349436)  0.8 
#>  4     1    61 (4654501 3220689, 4779626 3219592)  0.77
#>  5     1    88 (4654501 3220689, 4716357 3109639)  0.82
#>  6     2    15 (4590239 3111565, 4529015 3222722)  0.94
#>  7     2    69 (4590239 3111565, 4464101 3114115)  0.94
#>  8     2    87 (4590239 3111565, 4589159 3005137)  0.95
#>  9     2    88 (4590239 3111565, 4716357 3109639)  0.88
#> 10     3    22 (5486111 3937993, 5500145 3778167)  0.96
#> # ℹ 212 more rows
```
