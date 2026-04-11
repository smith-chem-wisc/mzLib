# Top-Down Testing Task Plan

This is the execution order for the top-down testing suite. It is written to get value quickly from MetaMorpheus results while leaving a clean slot for ProSightPD later.

## Phase 1: MetaMorpheus-Only Ingestion

1. Add a `Test/TopDownEngine/Parsers/MetaMorpheus/` test area.
2. Identify or create the MetaMorpheus result fixture files.
3. Write parser contract tests for the MetaMorpheus result shape.
4. Add malformed-input tests for missing columns, bad masses, and bad RT values.
5. Confirm the parsed records normalize into the shared top-down result model.

## Phase 2: Consensus Extraction

1. Define the unified `TopDownSearchResult` and `ConsensusFeature` expectations.
2. Build tiny synthetic records for mass and RT clustering.
3. Write tests for `k-of-n` support thresholds.
4. Write tests for ordering stability.
5. Write tests for reconciliation rules on mass, RT, and charge lists.

## Phase 3: Empirical Characterization

1. Pick one consensus feature fixture with clear signal boundaries.
2. Write tests for XIC width, isotopologue depth, and charge breadth.
3. Write tests for baseline noise floor and neighborhood density.
4. Write tests for output schema and summary statistics.

## Phase 4: Simulator MVP

1. Define the simulator input manifest format.
2. Build one fixed-seed synthetic dataset.
3. Write determinism tests for mzML and manifest output.
4. Write statistical sanity checks for envelope shape and noise.
5. Write provenance tests to ensure every peak is explainable.

## Phase 5: Engine Harness

1. Add one compact synthetic acceptance fixture.
2. Run the new engine against that fixture.
3. Assert recovery rate, mass accuracy, charge recovery, and quant accuracy.
4. Add a real-data smoke test using the MetaMorpheus consensus set.
5. Emit a metrics file for each run.

## Phase 6: ProSightPD Extension

1. Add parser tests for the SQL-backed ProSightPD output once the parser lands.
2. Add ProSightPD fixtures to the consensus set.
3. Re-run the consensus and harness tests with MetaMorpheus + ProSightPD together.
4. Tighten the acceptance thresholds against the larger consensus set.

## Done Means

- MetaMorpheus can be parsed and normalized end to end.
- The consensus extractor works on synthetic and MetaMorpheus-backed inputs.
- The simulator is deterministic under a fixed seed.
- The harness can report recovery and quant metrics on at least one real-data smoke set.
- ProSightPD has a clear integration path, even before its parser exists.
