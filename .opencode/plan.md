# Development Plan & Task Registry

## Sequence Converter Implementation Plan

This document tracks the implementation of a robust Sequence Converter system for mzLib that enables conversion of `IBioPolymerWithSetMods` sequences between different modification naming conventions (MetaMorpheus, UniProt, Unimod).

### High-Level Goals

1. **Improve ModificationConverter** - Leverage database references (RESID, PSI-MOD, Unimod IDs) for accurate modification matching instead of relying solely on name/mass heuristics
2. **Create SequenceConverter** - A service that converts full sequences between naming conventions  
3. **Support Use Cases**:
   - Retention time prediction (Chronologer) with standardized mod encoding
   - Database export with UniProt-compatible modifications for other search engines
   - Reading/writing results from various file formats with correct mod mapping
   - Getting full sequences with mass shifts for display/comparison

---

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

Format: `ID|TITLE|TASK_DOC|STATUS|DEPS`

Status values: `TODO`, `IN_PROGRESS`, `DONE`

<!-- TASKS:BEGIN -->
SC001|Build Modification Cross-Reference Index|tasks/SC001_mod_crossref_index.md|IN_PROGRESS|
SC002|Refactor ModificationConverter with Database IDs|tasks/SC002_refactor_mod_converter.md|TODO|SC001
SC003|Create SequenceConverter Core Class|tasks/SC003_sequence_converter_core.md|TODO|SC002
SC004|Add Mass Shift Sequence Formatting|tasks/SC004_mass_shift_formatting.md|TODO|SC003
SC005|Integrate with Chronologer|tasks/SC005_chronologer_integration.md|TODO|SC003
SC006|Integrate with ProteinDbWriter|tasks/SC006_proteindbwriter_integration.md|TODO|SC003
SC007|Add Unit Tests for SequenceConverter|tasks/SC007_unit_tests.md|TODO|SC003
SC008|Documentation and Examples|tasks/SC008_documentation.md|TODO|SC007
<!-- TASKS:END -->

---

## Architecture Overview

```
Omics/
+-- Modifications/
¦   +-- Modification.cs              # Existing - has DatabaseReference property
¦   +-- Mods.cs                      # Existing - static mod registries by convention
¦   +-- ModificationCrossRefIndex.cs # NEW - index for cross-referencing mods
¦   +-- SequenceConverter.cs         # NEW - main conversion service
¦
Readers/
+-- ExternalResults/SupportClasses/
¦   +-- ModificationConverter.cs     # REFACTOR - use cross-ref index
```

### Key Data Structures

The `Modification.DatabaseReference` property contains cross-reference IDs:
```
Dictionary<string, IList<string>> DatabaseReference:
  "RESID" -> ["AA0441"]
  "PSI-MOD" -> ["MOD:01624"]
  "Unimod" -> ["35"]
```

These can be used to link modifications across conventions:
- UniProt ptmlist.txt entries have `DR` (database reference) lines with RESID, PSI-MOD IDs
- Unimod modifications have `record_id` which becomes `Unimod: {id}` in DatabaseReference
- MetaMorpheus Mods.txt can have `DR` lines added

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
   dotnet test ... --filter "?..."
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
dotnet test C:/Users/Nic/Source/Repos/mzLib/mzLib/Test/Test.csproj --filter "FullyQualifiedName~Modification"
dotnet test C:/Users/Nic/Source/Repos/mzLib/mzLib/Test/Test.csproj --filter "FullyQualifiedName~SequenceConverter"

# Git status
git status
git log --oneline -10

# Commit
git add <files>
git commit -m "descriptive message"
