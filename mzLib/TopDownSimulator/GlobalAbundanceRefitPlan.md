# Global Abundance Refit Plan

## Goal

Improve agreement between simulated and real peak heights by adding a post-fit abundance refinement step that jointly adjusts proteoform abundances while keeping fitted shape parameters fixed.

## Current baseline

- `ParameterFitter` currently fits each proteoform independently from its extracted ground truth.
- `AbundanceFitter` solves a single-proteoform least-squares scale term with fixed `sigmaMz`, RT profile, and charge distribution.
- In crowded regions, independent fits can over/under-estimate abundance due to co-eluting signal shared across proteoforms.

## Proposed algorithm

1. Start from fitted `ProteoformModel[]` (initial abundances from `ParameterFitter`).
2. For each fitted proteoform, retain its corresponding `ProteoformGroundTruth` extraction window.
3. Build an iterative non-negative coordinate-descent refit:
   - For target proteoform `p`, compute denominator: `sum(phi_p^2)` over that truth window.
   - Compute numerator using residual with current contributions from all proteoforms at each sample point.
   - Update `A_p <- max(0, numerator / denominator)`.
4. Sweep all proteoforms for `N` iterations or until max relative abundance change is below tolerance.

This approximates NNLS behavior while remaining lightweight and compatible with existing forward-model primitives.

## API and code changes

- Add `TopDownSimulator/Fitting/GlobalAbundanceRefitter.cs`:
  - options record (`MaxIterations`, `ConvergenceTolerance`, `MinimumAbundance`, `Verbose`);
  - result record (`Models`, `IterationsCompleted`, `Converged`, `InitialResidualFraction`, `FinalResidualFraction`);
  - `Refit(IReadOnlyList<FittedProteoform> fitted, IReadOnlyList<ProteoformGroundTruth> truths, int minCharge, int maxCharge, double sigmaMz)`.
- Keep RT profile, charge distribution, monoisotopic mass, and identifier unchanged; only replace abundance.
- Integrate into `AnalysisExample.FitProteoforms(...)` via optional flag/env var.

## Validation strategy

- Unit tests:
  - recovery test on synthetic overlapping proteoforms where global refit improves abundance accuracy;
  - non-negativity test (all abundances remain `>= 0`);
  - improvement test that residual energy fraction decreases after refit.
- Workflow test harness:
  - log initial vs final residual fraction and iteration count for rep2 fract7 path.

## Performance considerations

- Use current `ForwardModel.Evaluate` for correctness-first implementation.
- Restrict each coordinate update to that proteoform's extracted truth window.
- Add early convergence stop.
- If needed later, cache basis evaluations in a second pass optimization.

## Rollout

- Default to enabled in `ExportRep2SliceAndQValueSimulations`.
- Allow opt-out via environment variable for A/B comparisons and quick rollback.
