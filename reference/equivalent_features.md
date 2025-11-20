# equivalent_features

Manually curated list of aminoacids and aminoacid-related metabolites
with corresponding metabolite identifiers (HMDB, KEGG, etc.)
irrespective of chirality.

## Usage

``` r
equivalent_features
```

## Format

An object of class `spec_tbl_df` (inherits from `tbl_df`, `tbl`,
`data.frame`) with 34 rows and 7 columns.

## Examples

``` r
data(equivalent_features)
head(equivalent_features)
#> # A tibble: 6 × 7
#>   TrivialName Class              OtherTrivialNames     chebi hmdb  pubchem INCHI
#>   <chr>       <chr>              <chr>                 <chr> <chr> <chr>   <chr>
#> 1 AABA        Aminoacids Related D-alpha-Aminobutyric… 2879… HMDB… 439691… InCh…
#> 2 Ac-Orn      Aminoacids Related N(5)-acetyl-L-ornith… 1654… HMDB… 439232… InCh…
#> 3 ADMA        Aminoacids Related Asymmetric dimethyla… 1792… HMDB… 123831… InCh…
#> 4 Ala         Aminoacids         D-Alanine, L-Alanine… 1557… HMDB… 71080,… InCh…
#> 5 alpha-AAA   Aminoacids Related Aminoadipic acid, D-… 3702… HMDB… 92136,… InCh…
#> 6 Arg         Aminoacids         D-Arginine, L-Argini… 1581… HMDB… 71070,… InCh…
```
