# Task: Refresh plan & task instructions

## Objective
Document the adjustments made on 2026-03-05 to keep the Sequence Converter implementation plan and associated task files aligned with the current code layout and workflow expectations.

## Work Completed
- Updated `.opencode/plan.md` startup routine to point contributors at `C:/Users/Nic/Source/Repos/mzLib` (the repo root).
- Clarified the architecture diagram so new files appear under `Omics/Modifications/Conversion/`.
- Adjusted the affected task specs (SC001, SC003, SC004, SC005, SC006) to reference the correct folders and integration points.

## Acceptance Criteria
- [x] Contributors landing on `.opencode/plan.md` see the correct working directory.
- [x] New tasks point at the actual `Conversion` folder for SequenceConverter components.
- [x] Chronologer and ProteinDbWriter task instructions reference existing classes/methods.
- [x] Changes logged in this checklist with status `DONE` in the plan registry.
