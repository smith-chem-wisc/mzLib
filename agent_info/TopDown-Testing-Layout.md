# Top-Down Testing Folder Layout

This is the proposed structure for the top-down testing workstream. It keeps the MetaMorpheus-first path unblocked while leaving a clear slot for ProSightPD once the SQL parser is available.

## Proposed Layout

```text
Test/TopDownEngine/
  Parsers/
    MetaMorpheus/
      MetaMorpheusResultParserTests.cs
      TestData/
    ProSightPD/
      ProSightPdParserTests.cs
      TestData/
  Consensus/
    ConsensusFeatureExtractorTests.cs
    TestData/
  Characterization/
    EmpiricalCharacterizationTests.cs
    TestData/
  Simulation/
    TopDownSimulatorTests.cs
    TestData/
  Harness/
    TopDownEngineAcceptanceTests.cs
    TestData/
  Regression/
    TopDownRegressionTests.cs
    TestData/
  Shared/
    TopDownTestFixtures.cs
    TopDownAssertions.cs
    TopDownTestPaths.cs
```

## Layout Rules

- Keep each stage isolated so parser failures do not mask consensus or simulator bugs.
- Put small gold fixtures next to the tests that consume them.
- Put shared helpers in `Shared/` only if two or more test areas need them.
- Keep ProSightPD tests present as a folder stub, but disabled or empty until the SQL parser exists.

## Fixture Conventions

- MetaMorpheus fixtures first: one small real result file plus one malformed file.
- Consensus fixtures: tiny synthetic records with hand-checkable mass/RT tolerances.
- Characterization fixtures: one feature-rich real-data sample and one minimal synthetic sample.
- Simulator fixtures: one canonical input manifest, one expected output manifest, and one fixed-seed mzML.
- Harness fixtures: a single compact dataset that can run in CI.

## Scope Boundary

- Phase 1 only consumes MetaMorpheus results.
- ProSightPD is a planned extension and should not block the first test harness.
- Once the SQL parser lands, add ProSightPD parser tests and then extend consensus coverage to include it.
