using System;
using System.Collections.Generic;
using System.Linq;
using Readers;

namespace TopDownSimulator.Extraction;

public sealed record MmResultRecord(
    string FileNameWithoutExtension,
    int PrecursorScanNumber,
    int Ms2ScanNumber,
    int PrecursorCharge,
    double MonoisotopicMass,
    double RetentionTime,
    double Score,
    string FullSequence,
    string? Accession,
    string Identifier);

/// <summary>
/// Loads MetaMorpheus `.psmtsv` search results into a compact record shape the
/// simulator can use for anchoring extraction and Phase 3 co-eluter discovery.
/// </summary>
public sealed class MmResultLoader
{
    public IReadOnlyList<MmResultRecord> Load(string psmTsvPath)
    {
        var file = new PsmFromTsvFile(psmTsvPath, new SpectrumMatchParsingParameters
        {
            ParseMatchedFragmentIons = false,
        });
        file.LoadResults();

        return file.Results
            .Where(p => p is not null)
            .Where(p => p.MonoisotopicMass > 0 && p.RetentionTime >= 0)
            .Select(p => new MmResultRecord(
                FileNameWithoutExtension: p.FileNameWithoutExtension,
                PrecursorScanNumber: p.PrecursorScanNum,
                Ms2ScanNumber: p.Ms2ScanNumber,
                PrecursorCharge: p.PrecursorCharge,
                MonoisotopicMass: p.MonoisotopicMass,
                RetentionTime: p.RetentionTime,
                Score: p.Score,
                FullSequence: p.FullSequence,
                Accession: p.Accession,
                Identifier: BuildIdentifier(p)))
            .OrderBy(p => p.FileNameWithoutExtension, StringComparer.OrdinalIgnoreCase)
            .ThenBy(p => p.RetentionTime)
            .ThenBy(p => p.MonoisotopicMass)
            .ToArray();
    }

    private static string BuildIdentifier(PsmFromTsv psm)
    {
        if (!string.IsNullOrWhiteSpace(psm.Accession))
            return $"{psm.Accession}:{psm.FullSequence}:{psm.Ms2ScanNumber}";

        return $"{psm.FileNameWithoutExtension}:{psm.Ms2ScanNumber}";
    }
}
