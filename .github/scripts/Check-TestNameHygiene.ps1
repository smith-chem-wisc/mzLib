<#
.SYNOPSIS
    Fails if any discovered test has a '.' in its name, which shreds the Visual Studio
    Test Explorer tree into phantom, un-runnable nodes.

.DESCRIPTION
    NUnit splices a test's name verbatim into its fully-qualified name
    (Namespace.Class.TestName). Visual Studio Test Explorer, vstest, and `dotnet test --filter`
    all build their hierarchy by splitting that FQN on '.'. A test named with dots in it is
    therefore chopped into intermediate tree nodes that look like namespaces but are not:

        Test.FileReadingTests.ProForma.ProFormaLayer1RoundTripTests.top-down_v2-4.2.1-03
        -> ...ProFormaLayer1RoundTripTests / "top-down_v2-4" / "2" / "1-03"

    Those middle nodes cannot be expanded, run, or filtered on, and the real test is effectively
    hidden. Dots reach test names through NUnit's explicit naming APIs -- TestCaseData.SetName(),
    [TestCase(TestName = "...")], and [SetName]-equivalent attributes -- typically when a name is
    interpolated from data (a spec section number, a version string, a file name).

    Fix: sanitize the name at the point it is built, e.g. id.Replace('.', '_'). Do not rename the
    underlying data; only the display name needs to be dot-free.

    Argument lists are exempt. Auto-generated names like ReadFile("oblm1.xml") or Scale(0.1d) carry
    dots inside the trailing "(...)", which the Test Explorer treats as a single leaf and does not
    split. Only the name to the left of the first '(' is checked.

.PARAMETER Project
    Path to the test .csproj to inspect.

.PARAMETER NoBuild
    Pass --no-build to dotnet test. Use in CI when a build step has already run.

.EXAMPLE
    pwsh .github/scripts/Check-TestNameHygiene.ps1 -Project mzLib/Test/Test.csproj -NoBuild
#>
[CmdletBinding()]
param(
    [string]$Project = "mzLib/Test/Test.csproj",
    [switch]$NoBuild
)

if (-not (Test-Path $Project)) {
    Write-Host "ERROR: test project not found: $Project" -ForegroundColor Red
    exit 1
}

Write-Host "Discovering tests in $Project ..."

$dotnetArgs = @('test', $Project, '--list-tests', '--verbosity', 'quiet')
if ($NoBuild) { $dotnetArgs += '--no-build' }

# Deliberately no `2>&1`: under PowerShell 7 redirecting a native command's stderr while
# $ErrorActionPreference is 'Stop' turns any stderr line into a terminating NativeCommandError,
# and `dotnet test` routinely writes warnings there. The test list goes to stdout regardless.
$output = & dotnet @dotnetArgs
$discoveryExit = $LASTEXITCODE

if ($discoveryExit -ne 0) {
    $output | Write-Host
    Write-Host "ERROR: test discovery failed (dotnet test --list-tests exited $discoveryExit)." -ForegroundColor Red
    exit 1
}

# `dotnet test --list-tests` prints a header, then one indented display name per test.
$names = $output |
    Where-Object { $_ -match '^\s{4,}\S' } |
    ForEach-Object { $_.Trim() }

if (@($names).Count -eq 0) {
    $output | Write-Host
    Write-Host "ERROR: no tests were discovered -- the output format may have changed." -ForegroundColor Red
    Write-Host "Refusing to pass vacuously." -ForegroundColor Red
    exit 1
}

# Strip the argument list; only the name proper is subject to the rule.
$offenders = @($names | Where-Object {
    $paren = $_.IndexOf('(')
    $head = if ($paren -ge 0) { $_.Substring(0, $paren) } else { $_ }
    $head.Contains('.')
} | Sort-Object -Unique)

Write-Host "Checked $(@($names).Count) discovered tests."

if ($offenders.Count -gt 0) {
    Write-Host ""
    Write-Host "ERROR: $($offenders.Count) test name(s) contain a '.' and will fragment the Test Explorer tree:" -ForegroundColor Red
    $offenders | ForEach-Object { Write-Host "  $_" -ForegroundColor Red }
    Write-Host ""
    Write-Host "These names come from SetName(...) or [TestCase(TestName = ...)]. Replace '.' with '_'"
    Write-Host "where the name is built (the test data itself does not need to change)."
    exit 1
}

Write-Host "OK: no test name contains a '.'." -ForegroundColor Green
exit 0
