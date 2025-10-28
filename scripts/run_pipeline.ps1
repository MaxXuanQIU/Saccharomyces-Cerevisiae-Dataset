<#
Runs the full pipeline:
1) Preprocess network + GO (build labels, select proteins)
2) Extract ESM-2 features and re-align
3) Inspect the resulting pickle
#>

$repoRoot = Split-Path $PSScriptRoot -Parent
Set-Location $repoRoot

if (Test-Path .\.venv\Scripts\Activate.ps1) {
  . .\.venv\Scripts\Activate.ps1
}

Write-Host "[1/3] Preprocessing network + GO..." -ForegroundColor Cyan
python -m src.data_processor
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host "[2/3] Extracting ESM-2 features..." -ForegroundColor Cyan
python -m src.esm_feature_extractor
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host "[3/3] Inspecting processed dataset..." -ForegroundColor Cyan
python -m src.data_inspection
exit $LASTEXITCODE
