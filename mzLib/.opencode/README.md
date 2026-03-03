# mzLib Agent Loop System

A simple, disciplined task orchestration system that maintains code quality through mandatory build/test/commit cycles.

## Core Principles

1. **One Task at a Time**: Focus on a single task until completion
2. **Build Before Commit**: Every commit must build successfully
3. **Test Before Done**: Tests must pass before marking complete
4. **Git Discipline**: All work is committed with descriptive messages
5. **Activity Tracking**: Every session logs what was accomplished

## Quick Start

### 1. Check Current Status
```powershell
# See what needs to be done
cat .opencode/plan.md

# See what was accomplished
cat .opencode/Activity.md

# Check git history
git log --oneline -10
```

### 2. Add a Task

Edit `.opencode/plan.md`, add between `<!-- TASKS:BEGIN -->` and `<!-- TASKS:END -->`:
```
M001|Fix Parsing Bug|tasks/001_fix_parser.md|TODO|
M002|Add New Feature|tasks/002_new_feature.md|TODO|M001
```

Create `.opencode/tasks/001_fix_parser.md`:
```markdown
# Task: Fix Parsing Bug

## Objective
Fix the issue in SpectrumMatchFromTsv parsing

## Implementation Steps
- [ ] Identify the parsing error
- [ ] Update parsing logic
- [ ] Add test case
- [ ] Verify fix

## Acceptance Criteria
- [ ] Build succeeds
- [ ] Tests pass
- [ ] Bug is fixed

## Verification
dotnet test Test/Test.csproj --filter "FullyQualifiedName~Readers"
```

### 3. Run the Agent Loop

```powershell
# Interactive mode (you do the work, loop verifies)
.\opencode\scripts\agent_loop.ps1

# Work on specific task
.\opencode\scripts\agent_loop.ps1 -TaskID M001

# Dry run to see what would happen
.\opencode\scripts\agent_loop.ps1 -DryRun
```

## How It Works

```
???????????????????????
? Read plan.md        ?
? Find next TODO task ?
???????????????????????
           ?
           ?
???????????????????????
? Load task file      ?
? Show instructions   ?
???????????????????????
           ?
           ?
???????????????????????
? Do the work         ?
? (you or AI)         ?
???????????????????????
           ?
           ?
???????????????????????
? Build solution      ?
? ? Fail ? Fix & retry?
? ? Pass ? Continue  ?
???????????????????????
           ?
           ?
???????????????????????
? Run tests           ?
? ? Fail ? Fix & retry?
? ? Pass ? Continue  ?
???????????????????????
           ?
           ?
???????????????????????
? Commit changes      ?
? (with task ID)      ?
???????????????????????
           ?
           ?
???????????????????????
? Mark task DONE      ?
? Log to Activity.md  ?
???????????????????????
```

## Configuration

Edit `.opencode/config.json`:

```json
{
  "agent_loop": {
    "auto_commit": true,              // Auto-commit on success
    "commit_on_success": true,        // Only commit if verified
    "build_before_commit": true,      // Always build first
    "test_before_commit": true,       // Always test first
    "require_tests_pass": true,       // Tests MUST pass
    "max_iterations_per_session": 50, // Safety limit
    "pause_on_error": true            // Stop on first error
  }
}
```

## Task File Format

```markdown
# Task: [Brief Title]

## Objective
What needs to be accomplished

## Implementation Steps
- [ ] Step 1: Do this
- [ ] Step 2: Do that
- [ ] Step 3: Verify it works

## Acceptance Criteria
- [ ] Builds successfully
- [ ] Tests pass
- [ ] Feature works as intended

## Verification Command
dotnet test ... --filter "..."

## Notes
Additional context or gotchas
```

## Session Workflow

### Start of Session
1. `cd C:/Users/Nic/Source/Repos/mzLib/mzLib/`
2. Review `Activity.md` for previous work
3. Check `git log --oneline -10`
4. Review `plan.md` for next task
5. Run agent loop

### During Session
- Work on one task at a time
- Build frequently
- Test your changes
- Commit working code

### End of Session
1. Ensure current work builds and tests pass
2. Commit any uncommitted work
3. Update `Activity.md` with progress
4. Update `plan.md` task status
5. **CRITICAL**: Only mark DONE if fully complete and committed

## Example Session

```powershell
PS> .\opencode\scripts\agent_loop.ps1

???????????????????????????????????????????????????????????
  mzLib Agent Loop Session
???????????????????????????????????????????????????????????
Working Directory: C:/Users/Nic/Source/Repos/mzLib/mzLib/
Max Iterations: 50
Auto Commit: True
Require Tests: True

?? Session Startup Checks:
   Current directory: C:/Users/Nic/Source/Repos/mzLib/mzLib
   Git status:
     M Readers/SpectrumMatchFromTsv.cs
   Recent commits:
     abc1234 M001: Update parsing methods
     def5678 M002: Add generic helpers

??? Iteration 1/50 ???
?? Working on: M003 - Fix GlycoPsmFromTsv parsing

?? Task has 3 steps:
   - Update constructor to use generic methods
   - Remove manual parsing code
   - Test with sample file

?? Please implement the task steps above.
   When done, press Enter to verify and commit...

[You do the work...]

?? Building solution...
? Build successful
?? Running tests...
? Tests passed
?? Committing changes...
? Committed: M003: Fix GlycoPsmFromTsv parsing

? Task M003 completed and verified!
```

## Git Integration

All commits follow this format:
```
MXXX: Brief description of changes

- Detailed bullet point 1
- Detailed bullet point 2
```

Example:
```
M001: Update SpectrumMatchFromTsv parsing methods

- Consolidated GetOptionalValue to handle nullable types
- Removed duplicate parsing code
- Updated GlycoPsmFromTsv and OsmFromTsv constructors
```

## Verification Commands

```powershell
# Build
dotnet build C:/Users/Nic/Source/Repos/mzLib/mzLib/mzLib.sln

# Test all
dotnet test C:/Users/Nic/Source/Repos/mzLib/mzLib/Test/Test.csproj

# Test specific
dotnet test Test/Test.csproj --filter "FullyQualifiedName~Readers"
dotnet test Test/Test.csproj --filter "FullyQualifiedName~Quantification"

# Git check
git status
git diff
git log --oneline -5
```

## Troubleshooting

### Build Fails
- Review build output
- Fix compilation errors
- Run loop again

### Tests Fail
- Check test output
- Fix failing tests
- Run loop again

### Task Won't Complete
- Ensure all steps are done
- Verify build and tests pass
- Check git status (changes committed?)

## Best Practices

1. **Small Tasks**: Break big work into small, testable chunks
2. **Test Often**: Run tests after each significant change
3. **Commit Often**: Small, focused commits are better
4. **Clear Messages**: Commit messages should explain "why"
5. **Update Activity**: Keep Activity.md current
6. **Don't Mark DONE Prematurely**: Only when truly complete

## Integration with OpenCode/Ollama

Your `.opencode/config.json` already configures Ollama models. The agent loop works with:

- **Manual Mode**: You implement, loop verifies
- **Ollama Mode**: Model implements, loop verifies
- **Claude Mode**: API implements, loop verifies

All modes follow the same verification rules:
- ? Build must pass
- ? Tests must pass
- ? Changes must commit
- ? Then mark DONE

---

**Remember**: The goal is maintainable, verified progress. Every task marked DONE should be production-ready code that builds, tests, and is committed.
