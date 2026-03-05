---
name: write-tests
description: Write NUnit tests for the mzLib codebase
---

## What I do
- Identify what to test from the target class or method
- Place the test file in the correct location under `mzLib/Test/`
- Write well-structured NUnit 4 tests following codebase conventions

## When to use me
Use this when you want tests written for a class, method, or feature.
Tell me the class or file to test, or I will infer it from context.

## Project conventions

**Framework & runner:** NUnit 4 (`NUnit` 4.1.0) + `NUnit3TestAdapter` + `Microsoft.NET.Test.Sdk`
**Target:** `net8.0-windows`, platform `x64`
**Test project:** `mzLib/Test/Test.csproj`
**Test Cases** Use `[Test]` for individual test methods, `[TestCase]` for parameterized tests, and `[TestFixture]` for test classes.`

### File placement
| Code under test | Test file location |
|---|---|
| `mzLib/Omics/Foo.cs` | `mzLib/Test/Omics/FooTests.cs` |
| `mzLib/Chromatography/Bar/Baz.cs` | `mzLib/Test/RetentionTimePrediction/BazTests.cs` |
| `mzLib/Readers/Thing.cs` | `mzLib/Test/FileReadingTests/ThingTests.cs` |
| Other single-library classes | `mzLib/Test/TestFoo.cs` (flat, legacy style) |

### Namespace
Match the test folder: `namespace Test.Omics`, `namespace Test.RetentionTimePrediction`, `namespace Test` (flat files).

### Class skeleton
```csharp
using NUnit.Framework;
using <ProductionNamespace>;
// ... other usings

namespace Test.<SubFolder>;

[TestFixture]
[System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
public class FooTests
{
    // Optional: per-test setup
    [SetUp]
    public void SetUp() { ... }

    // Optional: per-test teardown (dispose IDisposable objects)
    [TearDown]
    public void TearDown() { ... }

    [Test]
    public void MethodName_Scenario_ExpectedOutcome()
    {
        // Arrange
        ...

        // Act
        var result = ...;

        // Assert
        Assert.That(result, Is.EqualTo(expected));
    }
}
```

### Assertion style
Always use **constraint-model** assertions (`Assert.That`), **not** classic `Assert.AreEqual`:
```csharp
Assert.That(value, Is.EqualTo(42));
Assert.That(str,   Does.StartWith("-"));
Assert.That(str,   Does.Contain("PEPTIDE"));
Assert.That(obj,   Is.Not.Null);
Assert.That(coll,  Is.Empty);
Assert.That(d,     Is.GreaterThan(0));
Assert.That(d,     Is.EqualTo(expected).Within(0.01));  // floating-point tolerance
Assert.Throws<ArgumentNullException>(() => sut.Method(null));
```
> Exception: older flat test files use `NUnit.Framework.Legacy.ClassicAssert` — match the style of the file you're editing.

### Test method naming
`MethodName_Scenario_ExpectedOutcome` — imperative, describes the input condition and the expected result.

Examples:
- `ConvertFullSequence_UnmodifiedPeptide_ReturnsBaseSequence`
- `ChronologerFormatting_OxidizedMethionine_ConvertsToLowercase`
- `Enum_CanBeConvertedToString`

### Peptide construction (domain-specific)
```csharp
// Unmodified
var peptide = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>());

// With modifications — must add each mod to the dictionary
var mods = new Dictionary<string, Modification>
{
    { "Oxidation on M", ModificationConverter.AllModsKnown["Oxidation on M"] }
};
var peptide = new PeptideWithSetModifications("PEPTM[Oxidation on M]IDE", mods);

// Using Mods static registry (for conversion tests)
var peptide = new PeptideWithSetModifications(
    "AC[Common Fixed:Carbamidomethyl on C]DE",
    Mods.AllKnownProteinModsDictionary);
```

### IDisposable predictors
Wrap predictors in `using` or dispose in `[TearDown]`:
```csharp
using var predictor = new ChronologerRetentionTimePredictor();
```

### Region grouping (for larger test classes)
Group related tests with `#region` blocks:
```csharp
#region Happy-path Tests
...
#endregion

#region Edge Cases
...
#endregion

#region Error Handling
...
#endregion
```

### What to cover
For each class or method, write tests for:
1. **Happy path** — expected input produces expected output
2. **Edge cases** — empty/null inputs, boundary values, single-element collections
3. **Error handling** — exceptions thrown for invalid input, graceful failure modes
4. **Consistency** — same input always produces same output
5. **Integration** — components interacting (only if unit tests are insufficient)

## Process
1. Read the target source file to understand the public API
2. Identify the correct test file path and namespace
3. Check if a test file already exists — extend it rather than create a new one
4. Draft tests covering the scenarios above
5. Write the file (or edit the existing one)
6. Do NOT run the tests — leave that to the user
