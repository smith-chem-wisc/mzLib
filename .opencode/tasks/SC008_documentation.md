# Task: Documentation and Examples

## Objective
Create comprehensive documentation for the SequenceConverter feature, including usage examples, API reference, and integration guides.

## Implementation Steps

- [ ] Add XML documentation to all public classes and methods
- [ ] Create `docs/SequenceConverter.md` with usage guide
- [ ] Add code examples for common scenarios
- [ ] Document modification naming conventions
- [ ] Create migration guide for existing code
- [ ] Update README if appropriate
- [ ] Add inline examples in test files as documentation

## Documentation Files

### 1. API Reference (XML Comments)

All public methods should have comprehensive XML documentation:

```csharp
/// <summary>
/// Converts a full sequence from one modification naming convention to another.
/// </summary>
/// <param name="fullSequence">The full sequence including modification annotations.</param>
/// <param name="sourceConvention">The naming convention of the input sequence.</param>
/// <param name="targetConvention">The desired naming convention for the output.</param>
/// <param name="failureMode">How to handle modifications that cannot be converted.</param>
/// <returns>The sequence with modifications converted to the target convention.</returns>
/// <example>
/// <code>
/// var converter = SequenceConverter.Default;
/// string result = converter.ConvertFullSequence(
///     "PEPT[Common Fixed:Carbamidomethyl on C]IDE",
///     ModificationNamingConvention.MetaMorpheus,
///     ModificationNamingConvention.UniProt);
/// // Result: "PEPT[UniProt:Carboxymethyl cysteine on C]IDE"
/// </code>
/// </example>
public string ConvertFullSequence(...) { }
```

### 2. Usage Guide (Markdown)

Create `.opencode/docs/SequenceConverter.md`:

```markdown
# SequenceConverter Usage Guide

## Overview
The SequenceConverter allows conversion of peptide/oligo sequences between different 
modification naming conventions (MetaMorpheus, UniProt, Unimod).

## Quick Start
\`\`\`csharp
// Convert a sequence to UniProt format
var converter = SequenceConverter.Default;
string uniprotSeq = converter.ConvertFullSequence(
    peptide.FullSequence,
    ModificationNamingConvention.MetaMorpheus,
    ModificationNamingConvention.UniProt);
\`\`\`

## Common Use Cases

### 1. Export for Other Search Engines
...

### 2. Standardize Input from Multiple Sources
...

### 3. Display with Mass Shifts
...
```

### 3. Integration Examples

Document how to use with:
- Chronologer predictions
- ProteinDbWriter exports
- PSM file readers
- External tool input/output

## Content Sections

1. **Overview**: What the SequenceConverter does
2. **Quick Start**: 5-line example to get started
3. **Naming Conventions**: Explain MetaMorpheus, UniProt, Unimod formats
4. **API Reference**: All public methods with parameters
5. **Common Use Cases**: Step-by-step guides for typical scenarios
6. **Troubleshooting**: Common issues and solutions
7. **Performance Tips**: Best practices for large-scale conversions

## Acceptance Criteria
- [ ] All public methods have XML documentation
- [ ] Usage guide covers main scenarios
- [ ] Code examples compile and work
- [ ] Documentation is clear and accurate
- [ ] Changes committed to git

## Verification
- Review generated API docs
- Run code examples from documentation
- Have another developer review for clarity
