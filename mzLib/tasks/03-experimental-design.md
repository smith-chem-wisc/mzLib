# Task 3: Build Hardcoded SpikeIn-5mix-MS3 ExperimentalDesign

## Objective
Create a class implementing `IExperimentalDesign` that encodes the full channel-to-condition mapping for the SpikeIn-5mix-MS3 TMT dataset (14 files, 10 channels each = 140 IsobaricQuantSampleInfo objects).

## Background

### IExperimentalDesign interface (`MassSpectrometry/ExperimentalDesign/IExperimentalDesign.cs`)
```csharp
public interface IExperimentalDesign
{
    Dictionary<string, ISampleInfo[]> FileNameSampleInfoDictionary { get; }
}
```
Keys are file names (with extension). Values are arrays of ISampleInfo, one per channel. The array ORDER must match the order of values in `ISpectralMatch.QuantValues`.

### IsobaricQuantSampleInfo constructor (`MassSpectrometry/ExperimentalDesign/IsobaricQuantSampleInfo.cs`)
```csharp
public IsobaricQuantSampleInfo(
    string fullFilePathWithExtension,   // file path
    string condition,                    // e.g., "0.125", "0.5", "0.667", "1.0", "Reference"
    int biologicalReplicate,             // 1 or 2 (each condition appears twice per mixture)
    int technicalReplicate,              // 1, 2, or 3 (run within mixture)
    int fraction,                        // 0 (no fractionation)
    int plexId,                          // mixture number (1-5)
    string channelLabel,                 // "126", "127N", etc.
    double reporterIonMz,                // e.g., 126.12776
    bool isReferenceChannel)             // true for 126 and 131
```

## Channel-to-Condition Mapping Table

From the paper's Supplementary Figure S1.2. Each value is the relative UPS1 concentration. "Ref" = reference channel.

**TMT10-plex Channel Order**: 126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131

| Channel | Mixture 1 & 5 | Mixture 2 | Mixture 3 | Mixture 4 |
|---------|---------------|-----------|-----------|-----------|
| 126     | Ref           | Ref       | Ref       | Ref       |
| 127N    | 0.667         | 0.5       | 0.125     | 1         |
| 127C    | 0.125         | 1         | 0.667     | 0.5       |
| 128N    | 0.5           | 0.667     | 1         | 0.125     |
| 128C    | 1             | 0.125     | 0.5       | 0.667     |
| 129N    | 0.125         | 1         | 0.5       | 0.667     |
| 129C    | 0.5           | 0.667     | 0.125     | 1         |
| 130N    | 1             | 0.125     | 0.667     | 0.5       |
| 130C    | 0.667         | 0.5       | 1         | 0.125     |
| 131     | Ref           | Ref       | Ref       | Ref       |

**Biological replicate assignment**: Within each mixture, each non-reference condition appears exactly twice. The first occurrence (scanning 127N→130C) gets BiologicalReplicate=1, the second gets BiologicalReplicate=2.

For Mixture 1 & 5:
- 0.667: channels 127N (biorep 1), 130C (biorep 2)
- 0.125: channels 127C (biorep 1), 129N (biorep 2)
- 0.5: channels 128N (biorep 1), 129C (biorep 2)
- 1: channels 128C (biorep 1), 130N (biorep 2)

For Mixture 2:
- 0.5: channels 127N (biorep 1), 130C (biorep 2)
- 1: channels 127C (biorep 1), 129N (biorep 2)
- 0.667: channels 128N (biorep 1), 129C (biorep 2)
- 0.125: channels 128C (biorep 1), 130N (biorep 2)

For Mixture 3:
- 0.125: channels 127N (biorep 1), 129C (biorep 2)
- 0.667: channels 127C (biorep 1), 130N (biorep 2)
- 1: channels 128N (biorep 1), 130C (biorep 2)
- 0.5: channels 128C (biorep 1), 129N (biorep 2)

For Mixture 4:
- 1: channels 127N (biorep 1), 129C (biorep 2)
- 0.5: channels 127C (biorep 1), 130N (biorep 2)
- 0.125: channels 128N (biorep 1), 130C (biorep 2)
- 0.667: channels 128C (biorep 1), 129N (biorep 2)

## File Names (14 files)

```
161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_01
161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_02
161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture1_03
161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture2_01
161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture2_02
161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture2_03
161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture3_01
161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture3_02
161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture3_03
161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture4_02
161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture4_03
161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture5_01
161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture5_02
161117_SILAC_HeLa_UPS1_TMT10_SPS_MS3_Mixture5_03
```

**Note**: Mixture4_01 is missing. The file names above do NOT have extensions. The dictionary keys need to match what `Path.GetFileName(fullFilePath)` returns from the ISpectralMatch objects. Since SpectrumMatchFromTsv strips extensions, use keys WITHOUT extensions.

**Technical replicate assignment**: Within each mixture:
- `_01` = TechnicalReplicate 1
- `_02` = TechnicalReplicate 2
- `_03` = TechnicalReplicate 3

## TMT10-plex Reporter Ion m/z Values

```
126   → 126.12776
127N  → 127.12476
127C  → 127.13108
128N  → 128.12811
128C  → 128.13443
129N  → 129.13147
129C  → 129.13779
130N  → 130.13482
130C  → 130.14114
131   → 131.13818
```

## File to Create

### `C:/Users/Alex/Source/Repos/mzLib/mzLib/Test/Quantification/TestHelpers/SpikeInExperimentalDesign.cs`

```csharp
using MassSpectrometry;
using MassSpectrometry.ExperimentalDesign;

namespace Test.Quantification.TestHelpers;

public class SpikeInExperimentalDesign : IExperimentalDesign
{
    public Dictionary<string, ISampleInfo[]> FileNameSampleInfoDictionary { get; }

    public SpikeInExperimentalDesign()
    {
        FileNameSampleInfoDictionary = BuildDesign();
    }

    private static Dictionary<string, ISampleInfo[]> BuildDesign()
    {
        var design = new Dictionary<string, ISampleInfo[]>();

        // Define channel labels and m/z values
        string[] channelLabels = { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131" };
        double[] channelMzs = { 126.12776, 127.12476, 127.13108, 128.12811, 128.13443, 129.13147, 129.13779, 130.13482, 130.14114, 131.13818 };

        // Define condition mappings per mixture type
        // Index corresponds to channel order (0=126, 1=127N, ..., 9=131)
        // Value is (condition, biologicalReplicate) or null for reference
        // ... build the mapping tables and loop over files ...

        return design;
    }
}
```

Implement by defining the condition arrays for each mixture type, then looping over files to create the IsobaricQuantSampleInfo objects. Use a helper method like:

```csharp
private static ISampleInfo[] CreateChannelInfos(
    string fileName, int mixtureNumber, int techRep,
    string[] channelLabels, double[] channelMzs,
    (string condition, int bioRep)[] channelConditions)
```

Where `channelConditions` has 10 entries, one per channel.

## Verification
1. Build compiles
2. Write a unit test in `TmtSpikeInTests.cs` (or a temporary test) that:
   - Instantiates `SpikeInExperimentalDesign`
   - Asserts `FileNameSampleInfoDictionary.Count == 14`
   - For each file, asserts 10 ISampleInfo objects
   - Checks a few specific channels: e.g., Mixture1_01 channel 126 is Reference, channel 128C has condition "1"
   - Verifies reference channels (126, 131) have `IsReferenceChannel == true`
