# build.ps1 -- compile a Fortran-90 directory with Digital Visual Fortran 90 (DF98)
# usage:  pwsh -File build.ps1 -SrcDir <dir> -Exe  <path.exe> [-Opt 4] [-Debug]
param(
    [Parameter(Mandatory=$true)][string]$SrcDir,
    [Parameter(Mandatory=$true)][string]$Exe,
    [int]$Opt = 4,
    [switch]$Dbg,
    [string]$ExtraFiles = ""
)
$vs = "C:\Program Files\Microsoft Visual Studio"
$env:PATH    = "$vs\DF98\BIN;$vs\VC98\Bin;" + $env:PATH
$env:LIB     = "$vs\DF98\LIB;$vs\VC98\LIB"
$env:INCLUDE = "$vs\DF98\INCLUDE;$vs\VC98\INCLUDE"

# module files first (module names -> file names), then the program last
$mods = @()
$progs = @()
$others = @()
foreach ($f in (Get-ChildItem -Path $SrcDir -Filter *.f90 | Sort-Object Name)) {
    $txt = Get-Content -Raw -Path $f.FullName
    if ($txt -match '(?im)^\s*(program)\s+\w+') { $progs += $f.FullName }
    elseif ($txt -match '(?im)^\s*module\s+\w+\s*$') { $mods += $f.FullName }
    else { $others += $f.FullName }
}
$order = @($mods) + @($others) + @($progs)
if ($ExtraFiles -ne "") { $order += ($ExtraFiles -split ';') }

$flags = @("/nologo", "/c", "/optimize:$Opt", "/free", "/real_size:64")
if ($Dbg) { $flags = @("/nologo", "/c", "/check:all", "/debug:full", "/free", "/Od", "/real_size:64") }
$obj = @()
Push-Location $SrcDir
foreach ($src in $order) {
    $o = & "$vs\DF98\BIN\F90.EXE" @flags $src 2>&1
    if ($o) { $o | Where-Object { $_ -match 'Error|error' } | ForEach-Object { Write-Output $_ } }
    $base = [System.IO.Path]::GetFileNameWithoutExtension($src)
    $obj += (Join-Path $SrcDir ($base + ".obj"))
}
Pop-Location
$link = & "$vs\VC98\Bin\LINK.EXE" /nologo /OUT:$Exe @obj 2>&1
if ($link) { $link | ForEach-Object { Write-Output $_ } }
if (Test-Path $Exe) { Write-Output ("BUILD OK -> " + $Exe) } else { Write-Output "BUILD FAILED" }

