# Uncoarsen a District Matrix

After a cores analysis or other form of coarsening, sometimes you need
to be at the original geography level to be comparable. This takes in a
coarsened matrix and uncoarsens it to the original level

## Usage

``` r
redist.uncoarsen(plans, group_index)
```

## Arguments

- plans:

  A coarsened matrix of plans.

- group_index:

  The index used to coarsen the shape.

## Value

matrix
