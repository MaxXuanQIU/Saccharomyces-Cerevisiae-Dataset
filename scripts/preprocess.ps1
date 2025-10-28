<#
Step 1: Preprocess network + GO to build labels and select proteins
#>

$repoRoot = Split-Path $PSScriptRoot -Parent
Set-Location $repoRoot

if (Test-Path .\.venv\Scripts\Activate.ps1) {
  . .\.venv\Scripts\Activate.ps1
}

python -m src.data_processor
exit $LASTEXITCODE
