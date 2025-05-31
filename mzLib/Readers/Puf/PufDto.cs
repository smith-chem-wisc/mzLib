using System.CodeDom;
using System.Runtime.CompilerServices;
using MassSpectrometry;
using MzLibUtil;

namespace Readers.Puf;

/// <summary>
/// An entire Puf data file. 
/// </summary>
public class PufDataSet
{
    internal string Version { get; set; } = default!;
    internal List<PufMsMsExperiment> Experiments { get; set; } = new();
}

public class PufMsMsExperiment
{
    internal string Id { get; set; } = default!;
    internal string Source { get; set; } = default!;
    internal string Comment { get; set; } = default!;
    internal PufDataScan InstrumentData { get; set; } = default!;
    internal List<PufAnalysis> Analyses { get; set; } = new();
}

/// <summary>
/// All relevant data about the scan. Basically a .mgf or .msAlign file.
/// </summary>
internal class PufDataScan
{
    internal string FragmentationMethod { get; set; } = default!;
    internal string IonType { get; set; } = default!;
    internal List<PufIntact> Intacts { get; set; } = new();
    internal List<PufFragment> Fragments { get; set; } = new();

    // Conversion operator to MsDataScan
    public static explicit operator NeutralMassSpectrum(PufDataScan scan)
    {
        // Use all fragments as peaks
        var mzs = scan.Fragments.Select(f => f.MzMonoisotopic ?? 0).ToArray();
        var intensities = scan.Fragments.Select(f => f.Intensity ?? 0).ToArray();

        // Remove zero m/z
        var peaks = mzs.Zip(intensities, (mz, inten) => new { mz, inten })
                       .Where(p => p.mz > 0)
                       .ToArray();
        if (peaks.Length == 0)
            throw new InvalidOperationException("No valid peaks found in the scan.");
        mzs = peaks.Select(p => p.mz).ToArray();
        intensities = peaks.Select(p => p.inten).ToArray();
        var charges = peaks.Select(p => 1).ToArray(); // 1 for all fragments for neutral mass spectrum

        return new NeutralMassSpectrum(mzs, intensities, charges, true);
    }
}

/// <summary>
/// Ms1 peak and its deconvoluted values
/// </summary>
internal class PufIntact
{
    internal string Id { get; set; } = default!;
    internal double? MzMonoisotopic { get; set; }
    internal double? MzAverage { get; set; }
    internal double? MassMonoisotopic { get; set; }
    internal double? MassAverage { get; set; }
    internal double? Intensity { get; set; }
}

/// <summary>
/// Ms2 peak and its deconvoluted values. 
/// </summary>
internal class PufFragment
{
    internal string Id { get; set; } = default!;
    internal double? MzMonoisotopic { get; set; }
    internal double? MzAverage { get; set; }
    internal double? MassMonoisotopic { get; set; }
    internal double? MassAverage { get; set; }
    internal double? Intensity { get; set; }
}

/// <summary>
/// A Puf analysis is a search against a database with specific parameters.
/// </summary>
internal class PufAnalysis
{
    internal string Id { get; set; } = default!;
    internal string Type { get; set; } = default!;
    internal PufSearchParameters SearchParameters { get; set; } = new();
    internal PufResults Results { get; set; } = new();
}

/// <summary>
/// Parameters for a Puf search against a database.
/// </summary>
internal class PufSearchParameters
{
    internal PufDatabase Database { get; set; } = default!;
    internal List<string> PtmList { get; set; } = new();
    internal List<string> FixedModificationList { get; set; } = new();
    internal List<string> TerminalModificationList { get; set; } = new();
    internal string IntactTolerance { get; set; } = default!;
    internal string IntactToleranceUnit { get; set; } = default!;
    internal string IntactMassType { get; set; } = default!;
    internal string FragmentTolerance { get; set; } = default!;
    internal string FragmentToleranceUnit { get; set; } = default!;
    internal string FragmentMassType { get; set; } = default!;
    internal bool? DeltaM { get; set; }
    internal bool? Multiplexing { get; set; }
    internal bool? ITraq { get; set; }
    internal bool? Disulfide { get; set; }
    internal PufHitCriteria HitCriteria { get; set; } = default!;
}

/// <summary>
/// Name of database used
/// </summary>
internal class PufDatabase
{
    internal string InternalName { get; set; } = default!;
    internal string Display { get; set; } = default!;
}

/// <summary>
/// Criteria for hits in a Puf search.
/// </summary>
internal class PufHitCriteria
{
    internal int? Max { get; set; }
    internal int? MinimumMatchesNum { get; set; }
    internal int? MinimumMatchesPercent { get; set; }
    internal string Score { get; set; } = default!;
}

/// <summary>
/// Results of a Puf search against a database.
/// </summary>
internal class PufResults
{
    internal List<PufIdentificationCollection> HitLists { get; set; } = new();
    internal string Legend { get; set; } = default!;
}

internal class PufIdentificationCollection
{
    internal string IntactId { get; set; } = default!;
    internal List<PufIdentification> Hits { get; set; } = new();
}

internal class PufIdentification
{
    internal string Id { get; set; } = default!;
    internal string MatchingGeneId { get; set; } = default!;
    internal string ProteinForm { get; set; } = default!;
    internal string Description { get; set; } = default!;
    internal bool? Signalp { get; set; }
    internal bool? Propep { get; set; }
    internal PufMatchingSequence MatchingSequence { get; set; } = default!;
    internal int? SequenceLength { get; set; }
    internal double? TheoreticalMass { get; set; }
    internal double? MassDifferenceDa { get; set; }
    internal double? MassDifferencePpm { get; set; }
    internal PufScore Score { get; set; } = default!;
    internal List<PufMatchingFragment> MatchingFragments { get; set; } = new();
}

internal class PufMatchingSequence
{
    internal string Resid { get; set; } = default!;
    internal string Format { get; set; } = default!;
    internal string Display { get; set; } = default!;
}

internal class PufScore
{
    internal string PScore { get; set; } = default!;
    internal string Expected { get; set; } = default!;
    internal string McluckeyScore { get; set; } = default!;
}

internal class PufMatchingFragment
{
    internal string Name { get; set; } = default!;
    internal double? TheoreticalMass { get; set; }
    internal List<PufFragmentMatch> Matches { get; set; } = new();
}

internal class PufFragmentMatch
{
    internal string Id { get; set; } = default!;
    internal double? MassDifferenceDa { get; set; }
    internal double? MassDifferencePpm { get; set; }
}
