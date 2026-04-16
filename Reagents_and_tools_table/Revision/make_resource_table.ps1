$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = Split-Path -Parent $scriptDir
$lockPath = Join-Path $scriptDir "renv.lock"
$pipelinePath = Join-Path $scriptDir "run_pipeline.sh"
$tsvPath = Join-Path $scriptDir "renv_resource_table.tsv"
$csvPath = Join-Path $scriptDir "renv_resource_table.csv"
$fallbackCsvPath = Join-Path $scriptDir "renv_resource_table_updated.csv"

if (-not (Test-Path $lockPath)) {
    $lockPath = Join-Path $repoRoot "renv.lock"
}

if (-not (Test-Path $pipelinePath)) {
    $pipelinePath = Join-Path $repoRoot "run_pipeline.sh"
}

$lock = Get-Content $lockPath -Raw | ConvertFrom-Json
$biocVersion = $lock.Bioconductor.Version
$rVersion = $lock.R.Version

$pipelineContent = Get-Content $pipelinePath -Raw
$pipelinePackages = [System.Collections.Generic.HashSet[string]]::new([System.StringComparer]::OrdinalIgnoreCase)
$metaProVizDependencies = [System.Collections.Generic.HashSet[string]]::new([System.StringComparer]::OrdinalIgnoreCase)

function Get-PackageNameFromDependencyString {
    param (
        [string]$Dependency
    )

    if ([string]::IsNullOrWhiteSpace($Dependency)) {
        return $null
    }

    $name = ($Dependency -replace '\s*\(.*$', '').Trim()
    if ($name -eq 'R') {
        return $null
    }

    return $name
}

$packageBlock = [regex]::Match($pipelineContent, 'packages\s*<-\s*c\((?<list>[\s\S]*?)\)', [System.Text.RegularExpressions.RegexOptions]::Singleline)
if ($packageBlock.Success) {
    foreach ($match in [regex]::Matches($packageBlock.Groups['list'].Value, '"([^"]+)"')) {
        [void]$pipelinePackages.Add($match.Groups[1].Value)
    }
}

if ($pipelineContent -match 'install_github\("([^"]+)"\)') {
    $githubRepo = $Matches[1]
    $githubPackageName = Split-Path $githubRepo -Leaf
    [void]$pipelinePackages.Add($githubPackageName)
}

$metaProVizPackage = $lock.Packages.MetaProViz
if ($null -ne $metaProVizPackage) {
    foreach ($field in @('Depends', 'Imports', 'Suggests')) {
        $dependencies = $metaProVizPackage.$field
        if ($null -eq $dependencies) {
            continue
        }

        foreach ($dependency in $dependencies) {
            $dependencyName = Get-PackageNameFromDependencyString -Dependency $dependency
            if ($null -ne $dependencyName) {
                [void]$metaProVizDependencies.Add($dependencyName)
            }
        }
    }
}

$out = [System.Collections.Generic.List[string]]::new()
$out.Add("Reagent/Resource`tReference or Source`tIdentifier or Catalog Number`tLoadedBy_run_pipeline.sh`tMetaProViz_Dependency")
$out.Add("R`thttps://www.r-project.org`tR $rVersion`tFALSE`tFALSE")

foreach ($prop in $lock.Packages.PSObject.Properties) {
    $pkg = $prop.Value
    $name = $pkg.Package
    $version = $pkg.Version
    $source = $pkg.Source
    $rawUrl = $pkg.URL
    $url = if ($rawUrl) { ($rawUrl -split ",")[0].Trim() } else { "" }
    $loadedByPipeline = if ($pipelinePackages.Contains($name)) { "TRUE" } else { "FALSE" }
    $metaProVizDependency = if ($metaProVizDependencies.Contains($name)) { "TRUE" } else { "FALSE" }

    if ($source -eq "Bioconductor") {
        $ref = if ($url) { $url } else { "https://bioconductor.org/packages/$name" }
        $id  = "Bioconductor $biocVersion; $name $version"
    } elseif ($source -eq "Repository") {
        $ref = if ($url) { $url } else { "https://cran.r-project.org/package=$name" }
        $id  = "CRAN; $name $version"
    } elseif ($source -eq "GitHub") {
        $user = $pkg.RemoteUsername
        $repo = $pkg.RemoteRepo
        $ref  = "https://github.com/$user/$repo"
        $id   = "GitHub; $user/$repo $version"
    } else {
        $ref = $url
        $id  = "$source; $name $version"
    }

    $out.Add("$name`t$ref`t$id`t$loadedByPipeline`t$metaProVizDependency")
}

$out | Out-File $tsvPath -Encoding UTF8

$csvRows = Import-Csv $tsvPath -Delimiter "`t"
$tempCsvPath = [System.IO.Path]::GetTempFileName()

try {
    $csvRows | Export-Csv $tempCsvPath -NoTypeInformation -Encoding UTF8

    try {
        Move-Item $tempCsvPath $csvPath -Force
        $writtenCsvPath = $csvPath
    } catch {
        Copy-Item $tempCsvPath $fallbackCsvPath -Force
        Remove-Item $tempCsvPath -Force
        $writtenCsvPath = $fallbackCsvPath
    }
} finally {
    if (Test-Path $tempCsvPath) {
        Remove-Item $tempCsvPath -Force
    }
}

Write-Host "Done. Total rows (incl. header): $($out.Count)"
Write-Host "Top-level packages from run_pipeline.sh: $($pipelinePackages.Count)"
Write-Host "Direct MetaProViz dependencies from renv.lock: $($metaProVizDependencies.Count)"
Write-Host "CSV written to: $writtenCsvPath"
