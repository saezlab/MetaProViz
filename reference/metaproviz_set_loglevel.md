# Sets the log level for the package logger

Sets the log level for the package logger

## Usage

``` r
metaproviz_set_loglevel(level, target = "logfile")
```

## Arguments

- level:

  Character or class `loglevel`. The desired log level.

- target:

  Character, either 'logfile' or 'console'

## Value

Returns `NULL`.

## Examples

``` r
metaproviz_set_loglevel(logger::FATAL, target = "console")
```
