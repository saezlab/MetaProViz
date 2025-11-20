# Save the current package configuration

Save the current package configuration

## Usage

``` r
metaproviz_save_config(path = NULL, title = "default", local = FALSE)
```

## Arguments

- path:

  Path to the config file. Directories and the file will be created if
  don't exist.

- title:

  Save the config under this title. One config file might contain
  multiple configurations, each identified by a title.

- local:

  Save into a config file in the current directory instead of a user
  level config file. When loading, the config in the current directory
  has priority over the user level config.

## Value

Returns `NULL`.

## Examples

``` r
# restore the defaults and write them to the default config file:
metaproviz_reset_config()
metaproviz_save_config()
```
