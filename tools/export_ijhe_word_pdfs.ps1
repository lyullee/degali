param(
    [Parameter(Mandatory = $true)]
    [string]$Docx,
    [Parameter(Mandatory = $true)]
    [string]$Pdf
)

$ErrorActionPreference = "Stop"
$docxPath = (Resolve-Path -LiteralPath $Docx).Path
$pdfPath = [System.IO.Path]::GetFullPath($Pdf)
$word = $null
$document = $null
try {
    $word = New-Object -ComObject Word.Application
    $word.Visible = $false
    $word.DisplayAlerts = 0
    $word.AutomationSecurity = 3
    $document = $word.Documents.Open($docxPath, $false, $true, $false)
    $document.ExportAsFixedFormat($pdfPath, 17)
    [pscustomobject]@{ docx = $docxPath; pdf = $pdfPath; pages = $document.ComputeStatistics(2) } |
        ConvertTo-Json -Compress
}
finally {
    if ($null -ne $document) {
        $document.Close(0)
        [void][System.Runtime.InteropServices.Marshal]::FinalReleaseComObject($document)
    }
    if ($null -ne $word) {
        $word.Quit()
        [void][System.Runtime.InteropServices.Marshal]::FinalReleaseComObject($word)
    }
    [GC]::Collect()
    [GC]::WaitForPendingFinalizers()
}
