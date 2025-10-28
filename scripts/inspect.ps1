<#
Step 3: Inspect processed pickle; prints keys, shapes, and sample preview
#>

$repoRoot = Split-Path $PSScriptRoot -Parent
Set-Location $repoRoot

if (Test-Path .\.venv\Scripts\Activate.ps1) {
  . .\.venv\Scripts\Activate.ps1
}

python -m src.data_inspection
exit $LASTEXITCODE
