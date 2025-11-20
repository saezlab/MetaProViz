# Load the package configuration from a config file

Load the package configuration from a config file

## Usage

``` r
metaproviz_load_config(path = NULL, title = "default", user = FALSE, ...)
```

## Arguments

- path:

  Path to the config file.

- title:

  Load the config under this title. One config file might contain
  multple configurations, each identified by a title. If the title is
  not available the first section of the config file will be used.

- user:

  Force to use the user level config even if a config file exists in the
  current directory. By default, the local config files have prioroty
  over the user level config.

- ...:

  Passed to `yaml.load_file`.

## Value

Invisibly returns the config as a list.

## Examples

``` r
metaproviz_load_config()
```
