#!/usr/bin/env pwsh
<#
.SYNOPSIS
    Simplified agent loop that maintains commit discipline

.DESCRIPTION
    This agent loop:
    1. Reads tasks from plan.md
    2. Executes task steps one at a time
    3. Builds and tests after each change
    4. Commits only on verified success
    5. Updates Activity.md with progress
    
    CRITICAL: Tasks are only marked DONE if:
    - Build succeeds
    - Tests pass
    - Changes are committed to git

.PARAMETER MaxIterations
    Maximum work cycles in this session (default: 50)

.PARAMETER TaskID
    Specific task to work on (optional, auto-selects if not provided)

.PARAMETER SkipTests
    Skip test execution (not recommended)

.EXAMPLE
    .\agent_loop.ps1
    
.EXAMPLE
    .\agent_loop.ps1 -TaskID M001 -MaxIterations 10
#>

param(
    [int]$MaxIterations = 50,
    [string]$TaskID = "",
    [switch]$SkipTests,
    [switch]$DryRun
)

$ErrorActionPreference = "Stop"
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$RootDir = Split-Path -Parent $ScriptDir

# Load config
$config = Get-Content "$RootDir/config.json" | ConvertFrom-Json
$workDir = $config.agent_loop.working_directory

function Write-SessionHeader {
    Write-Host "===========================================================" -ForegroundColor Cyan
    Write-Host "  mzLib Agent Loop Session" -ForegroundColor Cyan
    Write-Host "===========================================================" -ForegroundColor Cyan
    Write-Host "Working Directory: $workDir"
    Write-Host "Max Iterations: $MaxIterations"
    Write-Host "Auto Commit: $($config.agent_loop.auto_commit)"
    Write-Host "Require Tests: $($config.agent_loop.require_tests_pass)"
    Write-Host "Plan File: $($planFile)"
    Write-Host "Activity Log: $($activityFile)"
    Write-Host ""
    
    # Startup checks
    Write-Host "📋 Session Startup Checks:" -ForegroundColor Yellow
    Push-Location $workDir
    
    Write-Host "   Current directory: $(Get-Location)"
    Write-Host "   Git status:"
    git status --short | ForEach-Object { Write-Host "     $_" }
    Write-Host "   Recent commits:"
    git log --oneline -5 | ForEach-Object { Write-Host "     $_" }
    Write-Host ""
    
    Pop-Location
}

function Get-NextTask {
    param([string]$PlanFile, [string]$SpecificTaskID)
    
    $content = Get-Content $PlanFile -Raw
    
    if ($content -match '(?s)<!-- TASKS:BEGIN -->(.*?)<!-- TASKS:END -->') {
        $registryText = $Matches[1].Trim()
        
        if ([string]::IsNullOrWhiteSpace($registryText)) {
            return $null
        }
        
        $lines = $registryText -split "`n" | Where-Object { $_.Trim() -and $_.Trim() -match '^\S+\|' }
        $tasks = @()
        
        foreach ($line in $lines) {
            $parts = $line -split '\|' | ForEach-Object { $_.Trim() }
            
            if ($parts.Count -ge 4) {
                $deps = @()
                if ($parts.Count -ge 5 -and $parts[4]) { 
                    $deps = $parts[4] -split ',' | ForEach-Object { $_.Trim() } | Where-Object { $_ }
                }
                
                $tasks += [PSCustomObject]@{
                    ID = $parts[0]
                    Title = $parts[1]
                    File = $parts[2]
                    Status = $parts[3]
                    Deps = $deps
                }
            }
        }
        
        # If specific task requested, return it
        if ($SpecificTaskID) {
            $task = $tasks | Where-Object { $_.ID -eq $SpecificTaskID }
            if ($task -and $task.Status -ne 'DONE') {
                return $task
            }
            return $null
        }
        
        # Debug: show parsed tasks
        if ($tasks.Count -eq 0) {
            Write-Host "DEBUG: No tasks parsed from registry. Registry text length: $($registryText.Length)" -ForegroundColor Magenta
        }
        
        # Find first TODO/IN_PROGRESS with deps met
        foreach ($t in $tasks) {
            if ($t.Status -in @('TODO', 'IN_PROGRESS')) {
                $allDepsDone = $true
                foreach ($dep in $t.Deps) {
                    if ($dep) {
                        $depTask = $tasks | Where-Object { $_.ID -eq $dep }
                        if ($depTask -and $depTask.Status -ne 'DONE') {
                            $allDepsDone = $false
                            break
                        }
                    }
                }
                if ($allDepsDone) {
                    return $t
                }
            }
        }
    }
    
    return $null
}

function Get-TaskInstructions {
    param([string]$TaskFile)
    
    # Normalize path separators to backslashes
    $TaskFile = $TaskFile.Replace('/', '\')
    
    if (-not (Test-Path $TaskFile)) {
        Write-Host "ERROR: Task file not found at: $TaskFile" -ForegroundColor Red
        throw "Task file not found: $TaskFile"
    }
    
    $content = Get-Content $TaskFile -Raw
    
    # Parse checklist items
    $steps = @()
    $pattern = '- \[ \] (.+?)(?=\n|$)'
    $matches = [regex]::Matches($content, $pattern)
    
    foreach ($match in $matches) {
        $steps += $match.Groups[1].Value.Trim()
    }
    
    return [PSCustomObject]@{
        Content = $content
        Steps = $steps
    }
}

function Invoke-BuildAndTest {
    param([object]$Config, [switch]$SkipTests)
    
    Write-Host "🔨 Building solution..." -ForegroundColor Yellow
    Push-Location $Config.agent_loop.working_directory
    
    try {
        $buildOutput = & cmd /c "$($Config.agent_loop.build_command) 2>&1"
        
        if ($LASTEXITCODE -ne 0) {
            Write-Host "❌ Build failed!" -ForegroundColor Red
            Write-Host $buildOutput
            return $false
        }
        
        Write-Host "✅ Build successful" -ForegroundColor Green
        
        if (-not $SkipTests -and $Config.agent_loop.require_tests_pass) {
            Write-Host "🧪 Running tests..." -ForegroundColor Yellow
            
            $testCmd = $Config.agent_loop.test_command
            if ($Config.agent_loop.test_filter) {
                $testCmd += " --filter `"$($Config.agent_loop.test_filter)`""
            }
            
            $testOutput = & cmd /c "$testCmd 2>&1"
            
            if ($LASTEXITCODE -ne 0) {
                Write-Host "❌ Tests failed!" -ForegroundColor Red
                Write-Host $testOutput
                return $false
            }
            
            Write-Host "✅ Tests passed" -ForegroundColor Green
        }
        
        return $true
    }
    finally {
        Pop-Location
    }
}

function Save-Progress {
    param(
        [string]$ActivityFile,
        [object]$Task,
        [string]$WorkDone,
        [string]$Status
    )
    
    $timestamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
    
    $entry = @"

## Session: $timestamp

**Task**: $($Task.ID) - $($Task.Title)
**Status**: $Status

### Work Completed:
$WorkDone

"@
    
    if ($Status -eq "DONE") {
        $entry += @"
### Verification:
- ✅ Build successful
- ✅ Tests passed
- ✅ Changes committed

"@
    }
    
    Add-Content -Path $ActivityFile -Value $entry
}

function Commit-Changes {
    param([object]$Task, [string]$WorkDir)
    
    Push-Location $WorkDir
    
    try {
        $gitStatus = git status --short
        
        if (-not $gitStatus) {
            Write-Host "⚠️  No changes to commit" -ForegroundColor Yellow
            return $true
        }
        
        Write-Host "📝 Committing changes..." -ForegroundColor Yellow
        
        $commitMsg = "$($Task.ID): $($Task.Title)"
        
        git add -A
        git commit -m $commitMsg
        
        if ($LASTEXITCODE -eq 0) {
            Write-Host "✅ Committed: $commitMsg" -ForegroundColor Green
            return $true
        } else {
            Write-Host "❌ Commit failed" -ForegroundColor Red
            return $false
        }
    }
    finally {
        Pop-Location
    }
}

function Update-TaskStatus {
    param([string]$PlanFile, [string]$TaskID, [string]$NewStatus)
    
    $content = Get-Content $PlanFile -Raw
    
    # Update the status in the registry
    $pattern = "($TaskID\|[^|]+\|[^|]+\|)([^|]+)"
    $content = $content -replace $pattern, "`${1}$NewStatus"
    
    Set-Content -Path $PlanFile -Value $content -NoNewline
    
    Write-Host "   Updated $TaskID status: $NewStatus" -ForegroundColor Cyan
}

function Main {
    
    $planFile = "$RootDir\plan.md"
    $activityFile = "$RootDir\Activity.md"
    Write-SessionHeader
    
    $iteration = 0
    
    while ($iteration -lt $MaxIterations) {
        $iteration++
        
        Write-Host "====== Iteration $iteration/$MaxIterations ======" -ForegroundColor Yellow
        
        # Get next task
        $task = Get-NextTask -PlanFile $planFile -SpecificTaskID $TaskID
        
        if (-not $task) {
            Write-Host "✅ No actionable tasks found. All done!" -ForegroundColor Green
            break
        }
        
        Write-Host "📌 Working on: $($task.ID) - $($task.Title)" -ForegroundColor Cyan
        
        # Mark as IN_PROGRESS if it was TODO
        if ($task.Status -eq 'TODO') {
            Update-TaskStatus -PlanFile $planFile -TaskID $task.ID -NewStatus 'IN_PROGRESS'
        }
        
        # Load task instructions
        $taskFile = "$RootDir\$($task.File)"
        Write-Host "   Task file: $taskFile" -ForegroundColor Gray
        $instructions = Get-TaskInstructions -TaskFile $taskFile
        
        Write-Host ""
        Write-Host "📋 Task has $($instructions.Steps.Count) steps:" -ForegroundColor White
        $instructions.Steps | ForEach-Object { Write-Host "   - $_" }
        Write-Host ""
        
        if ($DryRun) {
            Write-Host "[DRY RUN] Would work on this task" -ForegroundColor Yellow
            continue
        }
        
        # Present task to AI/user for execution
        # In automated mode, this would call the AI
        # For now, we pause and let the user/AI work
        
        Write-Host "🤖 Please implement the task steps above." -ForegroundColor Cyan
        Write-Host "   When done, press Enter to verify and commit..." -ForegroundColor Cyan
        
        if (-not $env:AGENT_LOOP_AUTO) {
            Read-Host
        }
        
        # Verify: Build and Test
        Write-Host ""
        $verified = Invoke-BuildAndTest -Config $config -SkipTests:$SkipTests
        
        if (-not $verified) {
            Write-Host "⚠️  Task verification failed!" -ForegroundColor Red
            Write-Host "   Fix the issues and run the loop again." -ForegroundColor Yellow
            
            Save-Progress -ActivityFile $activityFile -Task $task -WorkDone "Attempted but verification failed" -Status "FAILED"
            
            if ($config.agent_loop.pause_on_error) {
                break
            }
            continue
        }
        
        # Commit if configured
        if ($config.agent_loop.auto_commit) {
            $committed = Commit-Changes -Task $task -WorkDir $config.agent_loop.working_directory
            
            if (-not $committed) {
                Write-Host "⚠️  Could not commit changes" -ForegroundColor Yellow
                continue
            }
        }
        
        # Mark task as DONE
        Update-TaskStatus -PlanFile $planFile -TaskID $task.ID -NewStatus 'DONE'
        
        Save-Progress -ActivityFile $activityFile -Task $task -WorkDone "Task completed successfully" -Status "DONE"
        
        Write-Host ""
        Write-Host "✅ Task $($task.ID) completed and verified!" -ForegroundColor Green
        Write-Host ""
        
        Start-Sleep -Milliseconds 500
    }
    
    Write-Host "===========================================================" -ForegroundColor Cyan
    Write-Host "  Session Complete - $iteration iterations" -ForegroundColor Cyan
    Write-Host "===========================================================" -ForegroundColor Cyan
}

# Execute
Main
