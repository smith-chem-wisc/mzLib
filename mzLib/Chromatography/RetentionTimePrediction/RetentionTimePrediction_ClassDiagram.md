# Retention Time Prediction Infrastructure - Class Diagram

## Architecture Overview

```mermaid
classDiagram
    %% Core Abstractions
    class IRetentionPredictable {
        <<interface>>
        +BaseSequence string
        +FullSequence string
        +MonoisotopicMass double
        +FullSequenceWithMassShifts string
    }

    class IRetentionTimePredictor {
        <<interface>>
        +PredictorName string
        +SeparationType SeparationType
        +PredictRetentionTime(peptide) double?
        +GetFormattedSequence(peptide) string?
    }

    class RetentionTimePredictor {
        <<abstract>>
        +ModHandlingMode IncompatibleModHandlingMode
        #MinSequenceLength int
        #MaxSequenceLength int
        +PredictRetentionTime(peptide) double?
        #PredictCore(peptide, sequence)* double?
        +GetFormattedSequence(peptide)* string?
        #ValidateBasicConstraints(peptide) bool
    }

    %% Concrete Predictors
    class SSRCalc3Predictor {
        -_calculator SSRCalc3
        +PredictorName "SSRCalc3"
        +SeparationType HPLC
        #MinSequenceLength 4
    }

    class ChronologerPredictor {
        -_model Chronologer
        +PredictorName "Chronologer"
        +SeparationType HPLC
        #MaxSequenceLength 50
        +Dispose()
    }

    class CZEPredictor {
        -_columnLengthMeters double
        -_voltsPerMeter double
        +PredictorName "CZE"
        +SeparationType CZE
        #MinSequenceLength 6
        +ExperimentalElectrophoreticMobility(time) double
    }

    %% Supporting Components
    class SSRCalc3 {
        +ScoreSequence(sequence) double
        -Electric(sequence) double
        -Helicity1(sequence) double
        -Helicity2(sequence) double
    }

    class Chronologer {
        -_model TorchSharp.Module
        +Predict(tensor) Tensor
        +Dispose()
    }

    class CZECalculations {
        <<static>>
        +PredictedElectrophoreticMobility(seq, mass)$ double
        -PredictedChargeCorrected(seq)$ double
    }

    %% Enums
    class IncompatibleModHandlingMode {
        <<enumeration>>
        RemoveIncompatibleMods
        UsePrimarySequence
        ThrowException
        ReturnNull
    }

    class SeparationType {
        <<enumeration>>
        HPLC
        CZE
    }

    %% Relationships
    IRetentionTimePredictor <|.. RetentionTimePredictor
    RetentionTimePredictor <|-- SSRCalc3Predictor
    RetentionTimePredictor <|-- ChronologerPredictor
    RetentionTimePredictor <|-- CZEPredictor
    
    IRetentionTimePredictor ..> IRetentionPredictable : predicts
    RetentionTimePredictor ..> IncompatibleModHandlingMode : uses
    RetentionTimePredictor ..> SeparationType : uses
    
    SSRCalc3Predictor *-- SSRCalc3
    ChronologerPredictor *-- Chronologer
    CZEPredictor ..> CZECalculations
```

## Component Summary

### Core Architecture

| Component | Type | Purpose |
|-----------|------|---------|
| **IRetentionPredictable** | Interface | Defines what can be predicted (implemented by `PeptideWithSetModifications`) |
| **IRetentionTimePredictor** | Interface | Contract for all RT predictors |
| **RetentionTimePredictor** | Abstract Base | Common validation, modification handling, and flow control |

### Predictors

| Predictor | Separation | Algorithm | Min/Max Length | Key Feature |
|-----------|------------|-----------|----------------|-------------|
| **SSRCalc3** | HPLC | Krokhin 2006 | 4+ residues | Uses base sequence only |
| **Chronologer** | HPLC | Deep Learning | ≤50 residues | Supports specific PTMs |
| **CZE** | CZE | Mobility 2017 | 6+ residues | Electrophoretic migration |

### Modification Handling

| Mode | Behavior |
|------|----------|
| **RemoveIncompatibleMods** | Strip unsupported modifications, predict with remaining |
| **UsePrimarySequence** | Ignore all modifications, use base sequence only |
| **ThrowException** | Fail with descriptive error |
| **ReturnNull** | Return null if incompatible modifications present |

## Prediction Flow

```mermaid
flowchart TD
    A[IRetentionPredictable] --> B{Validate Basic<br/>Constraints}
    B -->|Invalid| C[Return null or throw]
    B -->|Valid| D{Format Sequence}
    D --> E{Incompatible<br/>Modifications?}
    E -->|No| F[PredictCore]
    E -->|Yes| G{Handling Mode}
    G -->|Remove| H[Filter Mods]
    G -->|Primary| I[Use Base Sequence]
    G -->|Throw| J[Throw Exception]
    G -->|Null| K[Return null]
    H --> F
    I --> F
    F --> L[Return RT]
    
    style A fill:#e1f5ff
    style F fill:#ffe1f5
    style L fill:#e1ffe1
```

## Design Patterns

### Template Method
```
RetentionTimePredictor.PredictRetentionTime():
  1. ValidateBasicConstraints()
  2. GetFormattedSequence()      [abstract]
  3. PredictCore()                [abstract]
  4. Return result
```

### Strategy (Enum-based)
```
IncompatibleModHandlingMode determines behavior:
  - RemoveIncompatibleMods → Filter and predict
  - UsePrimarySequence     → Use base only
  - ThrowException         → Fail explicitly  
  - ReturnNull             → Silent failure
```

## Extension Guide

To add a new predictor:

```csharp
public class NewPredictor : RetentionTimePredictor
{
    // 1. Set properties
    public override string PredictorName => "NewPredictor";
    public override SeparationType SeparationType => SeparationType.HPLC;
    protected override int MinSequenceLength => 5;
    
    // 2. Implement formatting
    public override string? GetFormattedSequence(
        IRetentionPredictable peptide, 
        out RetentionTimeFailureReason? failureReason)
    {
        // Convert to predictor-specific format
    }
    
    // 3. Implement prediction
    protected override double? PredictCore(
        IRetentionPredictable peptide,
        string? formattedSequence)
    {
        // Your algorithm here
    }
}
```

## Key Dependencies

```
Chromatography (this project)
  ↓
  IRetentionPredictable
  ↓
Proteomics
  ↓
  PeptideWithSetModifications : IRetentionPredictable
```

## References

- **SSRCalc3**: Krokhin OV et al., *Anal. Chem.* 2006, 78(22):7785-95
- **Chronologer**: Searle Lab, harmonized peptide libraries
- **CZE**: Krokhin OV et al., *Anal. Chem.* 2017, 89(3):2000-08
