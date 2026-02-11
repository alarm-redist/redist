# Iowa County File

This data contains geographic and demographic information on the 99
counties of the state of Iowa.

## Usage

``` r
data("iowa")
```

## Format

sf tibble containing columns for useful data related to the
redistricting process

- `fips`:

  The FIPS code for the county.

- `cd_2010`:

  The 2010 congressional district assignments.

- `pop`:

  The total population of the precinct, according to the 2010 Census.

- `white`:

  The non-Hispanic white population of the precinct.

- `black`:

  The non-Hispanic Black population of the precinct.

- `hisp`:

  The Hispanic population (of any race) of the precinct.

- `vap`:

  The voting-age population of the precinct.

- `wvap`:

  The white voting-age population of the precinct.

- `bvap`:

  The Black voting-age population of the precinct.

- `hvap`:

  The Hispanic voting-age population of the precinct.

- `tot_08`:

  Number of total votes for president in the county in 2008.

- `dem_08`:

  Number of votes for Barack Obama in 2008.

- `rep_08`:

  Number of votes for John McCain in 2008.

- `region`:

  The 28E agency regions for counties.

- `geometry`:

  The sf geometry column containing the geographic information.

## Examples

``` r
data(iowa)
print(iowa)
#> Simple feature collection with 99 features and 15 fields
#> Geometry type: MULTIPOLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 4081849 ymin: 2879102 xmax: 5834228 ymax: 4024957
#> Projected CRS: NAD83(HARN) / Iowa North (ftUS)
#> First 10 features:
#>     fips       name cd_2010    pop  white black hisp    vap  wvap bvap hvap
#> 1  19001      Adair       3   7682   7507    11  101   5957  5860    5   53
#> 2  19003      Adams       3   4029   3922     8   37   3180  3109    6   22
#> 3  19005  Allamakee       1  14330  13325   109  757  11020 10430   82  425
#> 4  19007  Appanoose       2  12887  12470    55  181   9993  9745   40   99
#> 5  19009    Audubon       4   6119   6007     9   37   4780  4714    5   27
#> 6  19011     Benton       1  26076  25387    93  275  19430 19068   49  155
#> 7  19013 Black Hawk       1 131090 109968 11493 4907 102594 89541 7677 2865
#> 8  19015      Boone       4  26306  25194   202  505  20027 19448  103  260
#> 9  19017     Bremer       1  24276  23459   186  239  18763 18242  155  137
#> 10 19019   Buchanan       1  20958  20344    59  243  15282 14979   32  128
#>    tot_08 dem_08 rep_08    region                       geometry
#> 1    4053   1924   2060     South MULTIPOLYGON (((4592338 328...
#> 2    2206   1118   1046     South MULTIPOLYGON (((4528041 315...
#> 3    7059   3971   2965 Northeast MULTIPOLYGON (((5422507 401...
#> 4    6176   2970   3086     South MULTIPOLYGON (((5032545 306...
#> 5    3435   1739   1634 Northwest MULTIPOLYGON (((4487363 341...
#> 6   13712   7058   6447 Southeast MULTIPOLYGON (((5246216 357...
#> 7   64775  39184  24662 Northeast MULTIPOLYGON (((5175640 369...
#> 8   13929   7356   6293   Central MULTIPOLYGON (((4741174 354...
#> 9   12871   6940   5741 Northeast MULTIPOLYGON (((5174636 379...
#> 10  10338   6050   4139 Northeast MULTIPOLYGON (((5302846 370...
```
