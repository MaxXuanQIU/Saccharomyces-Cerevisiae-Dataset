<#
Step 2: Extract ESM-2 embeddings and re-align dataset
#>

$repoRoot = Split-Path $PSScriptRoot -Parent
Set-Location $repoRoot

if (Test-Path .\.venv\Scripts\Activate.ps1) {
  . .\.venv\Scripts\Activate.ps1
}

python -m src.esm_feature_extractor
exit $LASTEXITCODE
