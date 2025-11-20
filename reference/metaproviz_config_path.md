# Current config file path of MetaProViz

Current config file path of MetaProViz

## Usage

``` r
metaproviz_config_path(user = FALSE)
```

## Arguments

- user:

  Logical: prioritize the user level config even if a config in the
  current working directory is available.

## Value

Character: path to the config file.

## Examples

``` r
metaproviz_config_path()
#> [1] "~/.config/MetaProViz/metaproviz.yml"
```
