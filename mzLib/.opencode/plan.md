# Development Plan & Task Registry

This file tracks all development milestones and tasks for the mzLib project.

## Session Startup Routine (READ THIS FIRST)

1. Run `pwd` to confirm you are in `C:/Users/Nic/Source/Repos/mzLib/mzLib/`
2. Read `Activity.md` in this directory to see what previous sessions accomplished
3. Run `git log --oneline -10` to see recent commits
4. Check the task list below - find the first task with status `TODO` whose dependencies are all `DONE`
5. Read the task's instruction file in `tasks/`
6. Work on that task. Build and test after each meaningful change.
7. Commit completed work with descriptive messages.

## End-of-Session Protocol (DO THIS BEFORE YOUR CONTEXT FILLS UP)

When you are nearing the end of your context window or finishing a task:
1. **Build**: `dotnet build "C:/Users/Nic/Source/Repos/mzLib/mzLib/mzLib.sln"`
2. **Test**: Run appropriate tests based on the task
3. **Commit** working changes: `git add <specific files> && git commit -m "descriptive message"`
4. **Update Activity.md**: Append a dated entry describing what you accomplished
5. **Update this file**: Change task status from `TODO` to `DONE` for completed tasks
6. **CRITICAL**: Do NOT mark a task DONE unless:
   - Code builds successfully
   - Tests pass (if applicable)
   - You have verified the functionality works
   - Changes are committed to git

---

## Task Registry (Machine Readable)

Format: `ID|TITLE|FILE|STATUS|DEPS`

Status values: `TODO`, `IN_PROGRESS`, `DONE`

<!-- TASKS:BEGIN -->

<!-- TASKS:END -->

---

## How to Add New Tasks

1. Add a line to the registry above (between BEGIN/END markers):
   ```
   M001|Task Title|tasks/001_description.md|TODO|
   ```

2. Create the task file in `.opencode/tasks/001_description.md`:
   ```markdown
   # Task: Task Title
   
   ## Objective
   What needs to be accomplished
   
   ## Implementation Steps
   - [ ] Step 1
   - [ ] Step 2
   - [ ] Step 3
   
   ## Acceptance Criteria
   - [ ] Builds successfully
   - [ ] Tests pass
   - [ ] Changes committed
   
   ## Verification Command
   dotnet test ... --filter "..."
   ```

3. Run the agent loop or work manually

---

## Quick Commands

```powershell
# Build
dotnet build C:/Users/Nic/Source/Repos/mzLib/mzLib/mzLib.sln

# Test all
dotnet test C:/Users/Nic/Source/Repos/mzLib/mzLib/Test/Test.csproj

# Test specific
dotnet test C:/Users/Nic/Source/Repos/mzLib/mzLib/Test/Test.csproj --filter "FullyQualifiedName~Quantification"

# Git status
git status
git log --oneline -10

# Commit
git add <files>
git commit -m "descriptive message"
