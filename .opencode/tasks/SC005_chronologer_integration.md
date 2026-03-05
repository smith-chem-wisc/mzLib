# Task: Integrate SequenceConverter with Chronologer

## Objective
Integrate the `SequenceConverter` with the Chronologer retention time prediction system to ensure peptide sequences are properly encoded for the neural network.

## Background

Chronologer expects peptide sequences in a specific encoding format. The `ChronologerEncoder.cs` class handles encoding peptides for the neural network. We need to ensure that:

1. Modifications are properly encoded (Chronologer has a specific vocabulary)
2. Unknown modifications are handled gracefully
3. Sequences from different sources (MetaMorpheus, external files) work correctly

## Implementation Steps

- [ ] Review `Chromatography/RetentionTimePrediction/Chronologer/ChronologerEncoder.cs`
- [ ] Identify how modifications are currently encoded
- [ ] Add `SequenceConverter` integration for standardizing input sequences
- [ ] Handle modifications not in Chronologer's vocabulary
- [ ] Add option to convert unknown mods to mass-based approximations
- [ ] Update any prediction clients that use Chronologer
- [ ] Add tests for converted sequences

## File Locations
- `Chromatography/RetentionTimePrediction/Chronologer/ChronologerEncoder.cs`
- `PredictionClients/` - any relevant prediction client code

## Key Integration Points

```csharp
// In ChronologerEncoder or ChronologerPredictor
public class ChronologerPredictor
{
    private readonly SequenceConverter _sequenceConverter;
    
    public float Predict(IBioPolymerWithSetMods peptide)
    {
        // Standardize modifications to Chronologer-compatible format
        var standardizedMods = _sequenceConverter.ConvertModifications(
            peptide.AllModsOneIsNterminus,
            ModificationNamingConvention.MetaMorpheus, // or custom Chronologer convention
            ConversionFailureMode.UseMassShift);
        
        // Encode and predict
        var encoded = Encode(peptide.BaseSequence, standardizedMods);
        return _model.Predict(encoded);
    }
}
```

## Modification Vocabulary

Chronologer likely supports a limited set of modifications. We need to:
1. Document which modifications are supported
2. Map common MetaMorpheus/UniProt mods to Chronologer vocabulary
3. Handle unsupported mods gracefully (mass shift encoding)

## Acceptance Criteria
- [ ] Chronologer works with sequences from different sources
- [ ] Unknown modifications are handled without errors
- [ ] Prediction accuracy is maintained
- [ ] Code builds without errors
- [ ] Changes committed to git

## Verification Commands
```powershell
# Build
dotnet build C:/Users/Nic/Source/Repos/mzLib/mzLib/mzLib.sln

# Test Chronologer
dotnet test C:/Users/Nic/Source/Repos/mzLib/mzLib/Test/Test.csproj --filter "FullyQualifiedName~Chronologer"
dotnet test C:/Users/Nic/Source/Repos/mzLib/mzLib/Test/Test.csproj --filter "FullyQualifiedName~RetentionTime"
```

## Notes
- Review the Chronologer paper for expected modification handling
- The model weights may have been trained with specific mod encodings
- Mass-based approximation may reduce prediction accuracy for rare mods
