# EPSG Table

This data contains NAD83 (HARN) EPSG codes for every U.S. state. Since
`redist` uses projected geometries, it is often a good idea to use
projections tailored to a particular state, rather than, for example, a
Mercator projection. Use these codes along with
[`sf::st_transform()`](https://r-spatial.github.io/sf/reference/st_transform.html)
to project your shapefiles nicely.

## Usage

``` r
data("EPSG")
```

## Format

named list containing EPSG codes for each U.S. state. Codes are indexed
by state abbreviations.

## Examples

``` r
data(EPSG)
EPSG$WA # 2855
#> [1] 2855
```
