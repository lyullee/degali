param(
    [string]$BrowserPath = "C:\Program Files\Google\Chrome\Application\chrome.exe"
)

$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot
$svg = [System.Uri]::new((Resolve-Path (Join-Path $root "IJHE_GRAPHICAL_ABSTRACT.svg"))).AbsoluteUri
$html = [System.Uri]::new((Resolve-Path (Join-Path $root "tmp\ijhe-graphical-abstract-print.html"))).AbsoluteUri
$png = Join-Path $root "outputs\ijhe-main-figure-images-2026-10-09-clean\figure-1.png"
$portalPdf = Join-Path $root "outputs\ijhe-graphical-abstract-2026-10-09.pdf"
$figurePdf = Join-Path $root "outputs\ijhe-figure-pdfs-2026-10-09\figure-1.pdf"
$pngProfile = Join-Path $root "tmp\chrome-ga-profile"
$pdfProfile = Join-Path $root "tmp\chrome-ga-pdf-profile"

if (-not (Test-Path -LiteralPath $BrowserPath)) {
    throw "Chrome executable not found: $BrowserPath"
}
New-Item -ItemType Directory -Force -Path (Split-Path -Parent $png), (Split-Path -Parent $figurePdf), $pngProfile, $pdfProfile | Out-Null

$pngArgs = @(
    "--headless=new", "--disable-gpu", "--hide-scrollbars", "--no-first-run",
    "--no-default-browser-check", "--user-data-dir=`"$pngProfile`"", "--window-size=1400,760",
    "--screenshot=`"$png`"", "`"$svg`""
)
$pngProcess = Start-Process -FilePath $BrowserPath -ArgumentList $pngArgs -Wait -PassThru -WindowStyle Hidden
if ($pngProcess.ExitCode -ne 0 -or -not (Test-Path -LiteralPath $png)) {
    throw "Graphical-abstract PNG export failed with exit code $($pngProcess.ExitCode)"
}

$pdfArgs = @(
    "--headless=new", "--disable-gpu", "--no-pdf-header-footer", "--print-to-pdf-no-header",
    "--no-first-run", "--no-default-browser-check", "--user-data-dir=`"$pdfProfile`"",
    "--run-all-compositor-stages-before-draw", "--print-to-pdf=`"$portalPdf`"", "`"$html`""
)
$pdfProcess = Start-Process -FilePath $BrowserPath -ArgumentList $pdfArgs -Wait -PassThru -WindowStyle Hidden
if ($pdfProcess.ExitCode -ne 0 -or -not (Test-Path -LiteralPath $portalPdf)) {
    throw "Graphical-abstract PDF export failed with exit code $($pdfProcess.ExitCode)"
}

Copy-Item -LiteralPath $portalPdf -Destination $figurePdf -Force
[pscustomobject]@{
    png = $png
    portal_pdf = $portalPdf
    figure_pdf = $figurePdf
} | ConvertTo-Json -Compress
