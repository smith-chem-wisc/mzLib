using System;
using System.Collections.Generic;
using System.Linq;

namespace TopDownSimulator.Extraction;

/// <summary>
/// Phase 3 helper that finds MetaMorpheus IDs near an anchor in retention time and
/// monoisotopic-mass space.
/// </summary>
public sealed class CoeluterFinder
{
    public IReadOnlyList<MmResultRecord> FindCoeluters(
        IReadOnlyList<MmResultRecord> results,
        MmResultRecord anchor,
        double rtHalfWidth,
        double massHalfWidthDa,
        bool sameFileOnly = true,
        bool includeAnchor = false)
    {
        if (rtHalfWidth < 0)
            throw new ArgumentOutOfRangeException(nameof(rtHalfWidth));
        if (massHalfWidthDa < 0)
            throw new ArgumentOutOfRangeException(nameof(massHalfWidthDa));

        return results
            .Where(candidate => !sameFileOnly || string.Equals(candidate.FileNameWithoutExtension, anchor.FileNameWithoutExtension, StringComparison.OrdinalIgnoreCase))
            .Where(candidate => includeAnchor || !IsSameRecord(candidate, anchor))
            .Where(candidate => Math.Abs(candidate.RetentionTime - anchor.RetentionTime) <= rtHalfWidth)
            .Where(candidate => Math.Abs(candidate.MonoisotopicMass - anchor.MonoisotopicMass) <= massHalfWidthDa)
            .OrderBy(candidate => Math.Abs(candidate.RetentionTime - anchor.RetentionTime))
            .ThenBy(candidate => Math.Abs(candidate.MonoisotopicMass - anchor.MonoisotopicMass))
            .ThenBy(candidate => candidate.RetentionTime)
            .ToArray();
    }

    private static bool IsSameRecord(MmResultRecord a, MmResultRecord b)
    {
        return string.Equals(a.FileNameWithoutExtension, b.FileNameWithoutExtension, StringComparison.OrdinalIgnoreCase)
               && a.Ms2ScanNumber == b.Ms2ScanNumber
               && a.PrecursorScanNumber == b.PrecursorScanNumber;
    }
}
