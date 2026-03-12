using System.Linq;
using System.Text;


namespace Omics.SequenceConversion;

/// <summary>
/// Serializes a <see cref="CanonicalSequence"/> into Chronologer format strings.
/// 
/// Chronologer uses single-character codes for modifications, with special
/// tokens for N-terminus and C-terminus states. This serializer converts
/// modifications based on their mass shift or resolved UNIMOD ID.
/// 
/// Output format example: "-PEPmTIDE_" (oxidized methionine)
/// - First character: N-terminus state ('-' = free)
/// - Middle: sequence with lowercase letters for modified residues
/// - Last character: C-terminus state ('_')
/// </summary>
public class ChronologerSequenceSerializer : SequenceSerializerBase
{
    /// <summary>
    /// Singleton instance with default settings.
    /// </summary>
    public static ChronologerSequenceSerializer Instance { get; } = new();

    /// <inheritdoc />
    public override string FormatName => ChronologerSequenceFormatSchema.Instance.FormatName;

    /// <inheritdoc />
    public override SequenceFormatSchema Schema => ChronologerSequenceFormatSchema.Instance;

    public ChronologerSequenceSerializer(IModificationLookup? lookup = null)
        : base(lookup ?? GlobalModificationLookup.Instance)
    {
    }

    #region Modification Mappings

    /// <summary>
    /// Maps UNIMOD IDs to Chronologer codes per target residue.
    /// Key: (UnimodId, TargetResidue), Value: Chronologer character
    /// </summary>
    private static readonly Dictionary<(int UnimodId, char Residue), char> UnimodToChronologer = new()
    {
        // Oxidation (UNIMOD:35)
        { (35, 'M'), 'm' },

        // Carbamidomethyl (UNIMOD:4)
        { (4, 'C'), 'c' },

        // Phosphorylation (UNIMOD:21)
        { (21, 'S'), 's' },
        { (21, 'T'), 't' },
        { (21, 'Y'), 'y' },

        // Acetyl (UNIMOD:1)
        { (1, 'K'), 'a' },

        // Succinyl (UNIMOD:64)
        { (64, 'K'), 'b' },

        // Ubiquitination (Dicarbamidomethyl on K, UNIMOD:121)
        { (121, 'K'), 'u' },

        // Methyl (UNIMOD:34)
        { (34, 'K'), 'n' },
        { (34, 'R'), 'q' },

        // Dimethyl (UNIMOD:36)
        { (36, 'K'), 'o' },
        { (36, 'R'), 'r' },

        // Trimethyl (UNIMOD:37)
        { (37, 'K'), 'p' },

        // TMT0 (UNIMOD:121)
        { (739, 'K'), 'z' },

        // TMT10 (UNIMOD:121)
        { (737, 'K'), 'x' },
    };

    /// <summary>
    /// Maps mass shifts (rounded to 2 decimal places) per target residue to Chronologer codes.
    /// Used as fallback when UNIMOD ID is not available.
    /// </summary>
    private static readonly Dictionary<(double MassRounded, char Residue), char> MassToChronologer = new()
    {
        // Oxidation on M: +15.995
        { (15.99, 'M'), 'm' },

        // Carbamidomethyl on C: +57.021
        { (57.02, 'C'), 'c' },

        // Phosphorylation: +79.966
        { (79.97, 'S'), 's' },
        { (79.97, 'T'), 't' },
        { (79.97, 'Y'), 'y' },

        // Acetylation on K: +42.011
        { (42.01, 'K'), 'a' },

        // Succinylation on K: +100.0
        { (100.02, 'K'), 'b' },

        // Ubiquitination on K (Dicarbamidomethyl): +114.04
        { (114.04, 'K'), 'u' },

        // Methylation on K/R: +14.016
        { (14.02, 'K'), 'n' },
        { (14.02, 'R'), 'q' },

        // Dimethylation on K/R: +28.031
        { (28.03, 'K'), 'o' },
        { (28.03, 'R'), 'r' },

        // Trimethylation on K: +42.047
        { (42.05, 'K'), 'p' },

        // TMT 0 on K: +224.12
        { (224.12, 'K'), 'z' },

        // TMT 10 on K: +229.16
        { (229.16, 'K'), 'x' },
    };

    /// <summary>
    /// N-terminal modification mappings by UNIMOD ID.
    /// </summary>
    private static readonly Dictionary<int, char> NTermUnimodCodes = new()
    {
        { 1, ChronologerSequenceFormatSchema.NTermAcetyl },   // Acetyl
        { 23, ChronologerSequenceFormatSchema.NTermPyroGlu },  // Pyro Glutamate from E
        { 24, ChronologerSequenceFormatSchema.NTermPyroGlu },  // Pyro Glutamate from Q
        { 26 , ChronologerSequenceFormatSchema.NTermCyclizedCamCys },  // cyclized CAM-Cys token.
        { 739, ChronologerSequenceFormatSchema.NTermTMT0 }, // TMT 0
        { 737, ChronologerSequenceFormatSchema.NTermTMT10 },  // TMT 10
    };

    /// <summary>
    /// N-terminal modification mappings by mass shift.
    /// </summary>
    private static readonly Dictionary<double, char> NTermMassCodes = new()
    {
        { 42.01, ChronologerSequenceFormatSchema.NTermAcetyl },     // Acetyl
        { -18.01, ChronologerSequenceFormatSchema.NTermPyroGlu },  // Pyro Glutamate from E
        { -17.03, ChronologerSequenceFormatSchema.NTermPyroGlu },  // Pyro Glutamate from Q
        { 39.99, ChronologerSequenceFormatSchema.NTermCyclizedCamCys },  // cyclized CAM-Cys token.
        { 224.12, ChronologerSequenceFormatSchema.NTermTMT0 },  // TMT 0
        { 229.16, ChronologerSequenceFormatSchema.NTermTMT10 }, // TMT 10
    };

    #endregion

    /// <inheritdoc />
    public override bool CanSerialize(CanonicalSequence sequence)
    {
        // Check sequence length constraint
        if (sequence.Length > ChronologerSequenceFormatSchema.MaxSequenceLength)
            return false;

        // Check for valid amino acids
        foreach (char aa in sequence.BaseSequence)
        {
            if (!ChronologerSequenceFormatSchema.CanonicalAminoAcids.Contains(aa))
                return false;
        }

        return true;
    }

    /// <inheritdoc />
    public override string? Serialize(
        CanonicalSequence sequence,
        ConversionWarnings? warnings = null,
        SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
    {
        warnings ??= new ConversionWarnings();

        // Validate sequence length
        if (sequence.Length > ChronologerSequenceFormatSchema.MaxSequenceLength)
        {
            return HandleError(warnings, mode, ConversionFailureReason.InvalidSequence,
                $"Sequence length {sequence.Length} exceeds Chronologer maximum of {ChronologerSequenceFormatSchema.MaxSequenceLength}.");
        }

        // Validate base sequence contains only canonical amino acids
        foreach (char aa in sequence.BaseSequence)
        {
            if (!ChronologerSequenceFormatSchema.CanonicalAminoAcids.Contains(aa))
            {
                return HandleError(warnings, mode, ConversionFailureReason.InvalidSequence,
                    $"Invalid amino acid '{aa}' not supported by Chronologer.");
            }
        }

        return base.Serialize(sequence, warnings, mode);
    }

    /// <inheritdoc />
    public override bool ShouldResolveMod(CanonicalModification mod)
    {
        return !mod.UnimodId.HasValue || !mod.EffectiveMass.HasValue;
    }

    protected override string? SerializeInternal(CanonicalSequence sequence, ConversionWarnings warnings, SequenceConversionHandlingMode mode)
    {
        var sb = new StringBuilder();

        // Build residue modifications lookup
        var residueMods = sequence.ResidueModifications
            .Where(m => m.ResidueIndex.HasValue)
            .ToDictionary(m => m.ResidueIndex!.Value, m => m);

        // Determine N-terminus token
        char nTermToken = DetermineNTerminusToken(sequence, residueMods, warnings, mode, out bool nTermHandled);
        
        // If N-term handling failed and we need to return null
        if (nTermToken == '\0')
        {
            if (mode == SequenceConversionHandlingMode.UsePrimarySequence)
                return $"{ChronologerSequenceFormatSchema.FreeNTerminus}{sequence.BaseSequence}{ChronologerSequenceFormatSchema.CTerminus}";
            if (mode == SequenceConversionHandlingMode.ReturnNull)
                return null;
        }

        sb.Append(nTermToken);

        // Process each residue
        for (int i = 0; i < sequence.BaseSequence.Length; i++)
        {
            char residue = sequence.BaseSequence[i];

            if (residueMods.TryGetValue(i, out var mod))
            {
                // Try to convert modification to Chronologer code
                char? chronologerCode = GetChronologerCode(mod, residue, warnings, mode);

                if (chronologerCode.HasValue)
                {
                    sb.Append(chronologerCode.Value);
                }
                else
                {
                    // Handle incompatible modification
                    switch (mode)
                    {
                        case SequenceConversionHandlingMode.ThrowException:
                            throw new SequenceConversionException(
                                $"Cannot convert modification '{mod.OriginalRepresentation}' on {residue} to Chronologer format.",
                                ConversionFailureReason.IncompatibleModifications,
                                new[] { mod.ToString() });

                        case SequenceConversionHandlingMode.ReturnNull:
                            return null;

                        case SequenceConversionHandlingMode.UsePrimarySequence:
                            return $"{ChronologerSequenceFormatSchema.FreeNTerminus}{sequence.BaseSequence}{ChronologerSequenceFormatSchema.CTerminus}";

                        case SequenceConversionHandlingMode.RemoveIncompatibleElements:
                        default:
                            warnings.AddWarning($"Removing incompatible modification: {mod.OriginalRepresentation} on {residue}");
                            warnings.AddIncompatibleItem(mod.ToString());
                            sb.Append(residue); // Use unmodified residue
                            break;
                    }
                }
            }
            else
            {
                sb.Append(residue);
            }
        }

        // Add C-terminus token
        sb.Append(ChronologerSequenceFormatSchema.CTerminus);

        return sb.ToString();
    }

    protected override string? GetModificationString(CanonicalModification mod, ConversionWarnings warnings, SequenceConversionHandlingMode mode)
    {
        throw new NotSupportedException("Chronologer serialization does not use bracketed modification strings.");
    }

    /// <summary>
    /// Determines the N-terminus token based on N-terminal modifications.
    /// </summary>
    private char DetermineNTerminusToken(
        CanonicalSequence sequence,
        Dictionary<int, CanonicalModification> residueMods,
        ConversionWarnings warnings,
        SequenceConversionHandlingMode mode,
        out bool handled)
    {
        handled = true;
        var nTermMod = sequence.NTerminalModification;

        // Check for PyroGlu at position 0
        if (residueMods.TryGetValue(0, out var firstResidueMod))
        {
            char firstResidue = sequence.BaseSequence[0];
            
            // Check for PyroGlu from E/Q at N-terminus
            if (firstResidue == 'E' || firstResidue == 'Q')
            {
                var code = GetChronologerCode(firstResidueMod, firstResidue, warnings, mode);
                if (code == 'e')
                {
                    // PyroGlu detected - use appropriate N-term token
                    return firstResidue == 'E' 
                        ? ChronologerSequenceFormatSchema.NTermPyroGlu 
                        : ChronologerSequenceFormatSchema.NTermPyroGlu;
                }
            }
            
            // Check for cyclized CAM-Cys at position 0
            if (firstResidue == 'C')
            {
                var code = GetChronologerCode(firstResidueMod, firstResidue, warnings, mode);
                if (code == 'd')
                {
                    return ChronologerSequenceFormatSchema.NTermCyclizedCamCys;
                }
            }
        }

        // Check explicit N-terminal modification
        if (nTermMod.HasValue)
        {
            var mod = nTermMod.Value;

            // Try UNIMOD ID first
            if (mod.UnimodId.HasValue && NTermUnimodCodes.TryGetValue(mod.UnimodId.Value, out char unimodCode))
            {
                return unimodCode;
            }

            // Try mass-based lookup
            if (mod.EffectiveMass.HasValue)
            {
                double mass = Math.Round(mod.EffectiveMass.Value, 2);
                if (NTermMassCodes.TryGetValue(mass, out char massCode))
                {
                    return massCode;
                }
            }

            // N-terminal modification not recognized
            handled = false;
            warnings.AddIncompatibleItem($"N-terminal: {mod.OriginalRepresentation}");
            warnings.AddWarning($"Unrecognized N-terminal modification: {mod.OriginalRepresentation}");
            
            if (mode == SequenceConversionHandlingMode.ThrowException)
            {
                throw new SequenceConversionException(
                    $"Cannot convert N-terminal modification '{mod.OriginalRepresentation}' to Chronologer format.",
                    ConversionFailureReason.IncompatibleModifications,
                    new[] { mod.ToString() });
            }

            return mode == SequenceConversionHandlingMode.ReturnNull 
                ? '\0' 
                : ChronologerSequenceFormatSchema.FreeNTerminus;
        }

        // Return free N-terminus by default
        return ChronologerSequenceFormatSchema.FreeNTerminus;
    }

    /// <summary>
    /// Gets the Chronologer single-character code for a modification.
    /// </summary>
    private char? GetChronologerCode(CanonicalModification mod, char targetResidue, ConversionWarnings warnings, SequenceConversionHandlingMode mode)
    {
        // Try UNIMOD ID first (most reliable)
        if (mod.UnimodId.HasValue)
        {
            if (UnimodToChronologer.TryGetValue((mod.UnimodId.Value, targetResidue), out char unimodCode))
            {
                return unimodCode;
            }
        }

        // Try mass-based lookup (fallback)
        if (mod.EffectiveMass.HasValue)
        {
            double massRounded = Math.Round(mod.EffectiveMass.Value, 2);
            if (MassToChronologer.TryGetValue((massRounded, targetResidue), out char massCode))
            {
                return massCode;
            }
        }

        // Could not map this modification
        return null;
    }

}
