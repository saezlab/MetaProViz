# Restore the built-in default values of all config parameters of MetaProViz

Restore the built-in default values of all config parameters of
MetaProViz

## Usage

``` r
metaproviz_reset_config(save = NULL, reset_all = FALSE)
```

## Arguments

- save:

  If a path, the restored config will be also saved to this file. If
  TRUE, the config will be saved to the current default config path (see
  [`metaproviz_config_path`](metaproviz_config_path.md)).

- reset_all:

  Reset to their defaults also the options already set in the R options.

## Value

The config as a list.

## See also

[`metaproviz_load_config`](metaproviz_load_config.md)`, `[`metaproviz_save_config`](metaproviz_save_config.md)

## Examples

``` r
# restore the defaults and write them to the default config file:
metaproviz_reset_config()
metaproviz_save_config()
```
