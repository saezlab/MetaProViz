# Pipe operator

See
`magrittr::`[`%>%`](https://magrittr.tidyverse.org/reference/pipe.html)
for details.

## Usage

``` r
lhs %>% rhs
```

## Arguments

- lhs:

  A value or the magrittr placeholder.

- rhs:

  A function call using the magrittr semantics.

## Value

The result of calling `rhs(lhs)`.

## Examples

``` r
c(1, 4, 9) %>% sqrt()
#> [1] 1 2 3
```
