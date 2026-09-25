using System.Globalization;
using System.IO.Compression;
using System.Text.RegularExpressions;
using System.Xml;
using MzIdentML;
using MzLibUtil;
using Omics;
using Omics.Modifications;
using Omics.SequenceConversion;
using Readers.ExternalResults.IndividualResultRecords;

namespace Readers.ExternalResults.ResultFiles;

/// <summary>
/// An mzIdentML file (.mzid), or a gzip-compressed one (.mzid.gz) read without decompressing it to disk.
/// Every schema version MzidIdentifications reads is supported, and every SpectrumIdentificationList in the
/// file is read.
/// </summary>
/// <remarks>
/// Each SpectrumIdentificationItem becomes one <see cref="MzIdentMLRecord"/>, except the items a record cannot
/// represent. Those are listed in <see cref="SkippedMatches"/> rather than failing the file:
/// <list type="bullet">
/// <item>crosslink items (MS:1002511), whose two peptides are not one linear match;</item>
/// <item>items whose peptide carries a modification with no UNIMOD accession, or one the Unimod catalog does
/// not resolve at that position, since a partial modification set would make FullSequence wrong;</item>
/// <item>items whose peptide carries a substitution, or two modifications at one position;</item>
/// <item>items whose peptide reference does not resolve.</item>
/// </list>
/// </remarks>
public class MzIdentMLResultFile : ResultFile<MzIdentMLRecord>
{
    // MS:1002354 PSM-level q-value and its is_a children in the shipped psi-ms.obo
    private static readonly HashSet<string> QValueAccessions =
    [
        "MS:1002354", // PSM-level q-value
        "MS:1002054", // MS-GF:QValue
        "MS:1002723", // MSPathFinder:QValue
        "MS:1003125", // ProSight:spectral Q-value
    ];

    private const string CrosslinkItemAccession = "MS:1002511";
    private static readonly Regex NativeIdIndex = new(@"(^|\s)index=(\d+)($|\s)", RegexOptions.Compiled);
    private static readonly Regex NativeIdFile = new(@"(^|\s)file=(\S+)", RegexOptions.Compiled);

    public override SupportedFileType FileType => FilePath.ParseFileType();
    public override Software Software { get; set; }

    private List<MzIdentMLSkippedMatch> _skippedMatches = [];
    private string? _loadedFrom;

    /// <summary>
    /// The items that were not turned into records, and why. Reading this forces the lazy load.
    /// </summary>
    public IReadOnlyList<MzIdentMLSkippedMatch> SkippedMatches
    {
        get
        {
            EnsureLoaded();
            return _skippedMatches;
        }
    }

    public MzIdentMLResultFile(string filePath) : base(filePath, Software.MzIdentML) { }

    /// <summary>Constructor used to initialize from the factory method.</summary>
    public MzIdentMLResultFile() : base()
    {
        Software = Software.MzIdentML;
    }

    private void EnsureLoaded()
    {
        _ = Results;
    }

    /// <summary>
    /// Reads the file. A file is read once per <see cref="ResultFile{TResult}.FilePath"/>: ResultFile reloads
    /// whenever Results is empty, and a file whose every item is skipped (a crosslink search, say) would
    /// otherwise be parsed again on every access.
    /// </summary>
    /// <exception cref="MzLibException">The file is not a readable mzIdentML document.</exception>
    public override void LoadResults()
    {
        if (_loadedFrom == FilePath)
        {
            return;
        }

        MzidIdentifications identifications;
        using (var file = new FileStream(FilePath, FileMode.Open, FileAccess.Read, FileShare.Read))
        using (var stream = FilePath.EndsWith(".gz", StringComparison.OrdinalIgnoreCase)
                   ? new GZipStream(file, CompressionMode.Decompress)
                   : (Stream)file)
        {
            try
            {
                identifications = new MzidIdentifications(stream);
            }
            catch (Exception e) when (e is InvalidOperationException or XmlException or InvalidDataException)
            {
                throw new MzLibException($"Could not read mzIdentML file '{FilePath}': {e.Message}", e);
            }
        }

        var spectraData = new Dictionary<string, MzidSpectraData>();
        foreach (var entry in identifications.GetSpectraData())
        {
            if (entry.Id != null)
            {
                spectraData.TryAdd(entry.Id, entry);
            }
        }

        var lookup = new UnimodModificationLookup();
        var results = new List<MzIdentMLRecord>();
        var skipped = new List<MzIdentMLSkippedMatch>();

        foreach (var match in identifications.GetSpectrumMatches())
        {
            string? reason = ReasonToSkip(match, lookup, out var mods);
            if (reason != null)
            {
                skipped.Add(new MzIdentMLSkippedMatch(match.SpectrumIdentificationItemId, match.SpectrumId, reason));
                continue;
            }

            results.Add(ToRecord(match, mods, spectraData));
        }

        _skippedMatches = skipped;
        Results = results;
        _loadedFrom = FilePath;
    }

    /// <summary>
    /// mzIdentML cannot be regenerated from these records, so writing is not supported. MetaMorpheus writes
    /// mzIdentML from its own search results.
    /// </summary>
    /// <exception cref="NotSupportedException">Always.</exception>
    public override void WriteResults(string outputPath) =>
        throw new NotSupportedException("Writing mzIdentML is not supported; MzIdentMLResultFile only reads it.");

    private static MzIdentMLRecord ToRecord(MzidSpectrumMatch match, Dictionary<int, Modification> mods,
        Dictionary<string, MzidSpectraData> spectraData)
    {
        var spectraFile = match.SpectraData;
        var fileId = match.SpectrumId == null ? null : NativeIdFile.Match(match.SpectrumId);
        if (fileId is { Success: true } && spectraData.TryGetValue(fileId.Groups[2].Value, out var named))
        {
            spectraFile = named;
        }

        string location = spectraFile?.Location ?? spectraFile?.Name ?? string.Empty;
        string baseSequence = match.PeptideSequence!;

        return new MzIdentMLRecord
        {
            Match = match,
            OneBasedScanNumber = OneBasedScanNumberOf(match.SpectrumId),
            SpectrumId = match.SpectrumId ?? string.Empty,
            SpectrumTitle = match.SpectrumTitle,
            SpectraFileLocation = location,
            // locations are written on the submitter's machine, so either separator ends a directory on any OS
            FileNameWithoutExtension = PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(
                location[(location.LastIndexOfAny(['/', '\\']) + 1)..]),
            BaseSequence = baseSequence,
            AllModsOneIsNterminus = mods,
            FullSequence = IBioPolymerWithSetMods.DetermineFullSequence(baseSequence, mods),
            Accession = string.Join('|', match.PeptideEvidence
                .Select(e => e.DBSequenceAccession)
                .Where(a => !string.IsNullOrEmpty(a))
                .Distinct()),
            IsDecoy = match.PeptideEvidence.Count > 0 && match.PeptideEvidence.All(e => e.IsDecoy),
            ChargeState = match.ChargeState,
            ExperimentalMassToCharge = match.ExperimentalMassToCharge,
            CalculatedMassToCharge = match.CalculatedMassToCharge,
            Rank = match.Rank,
            PassThreshold = match.PassThreshold,
            Scores = ScoresOf(match),
            QValue = QValueOf(match),
        };
    }

    /// <summary>
    /// The nativeID's scan number: "scan=N" gives N; "index=N", a zero-based position in a peak list, gives
    /// N + 1, as CasanovoMzTabFile does; anything else gives -1.
    /// </summary>
    internal static int OneBasedScanNumberOf(string? spectrumId)
    {
        if (string.IsNullOrEmpty(spectrumId))
        {
            return -1;
        }

        var scan = Mzml.nativeIdScanNumberParser.Match(spectrumId);
        if (scan.Success && int.TryParse(scan.Groups[2].Value, NumberStyles.None, CultureInfo.InvariantCulture, out int scanNumber))
        {
            return scanNumber;
        }

        var index = NativeIdIndex.Match(spectrumId);
        if (index.Success && int.TryParse(index.Groups[2].Value, NumberStyles.None, CultureInfo.InvariantCulture, out int position))
        {
            return position + 1;
        }

        return -1;
    }

    /// <summary>
    /// Null when the item can be a record, with its modifications in <paramref name="mods"/> keyed the
    /// AllModsOneIsNterminus way; otherwise why it cannot.
    /// </summary>
    private static string? ReasonToSkip(MzidSpectrumMatch match, UnimodModificationLookup lookup,
        out Dictionary<int, Modification> mods)
    {
        mods = [];

        if (match.ItemCvParams.Any(cv => cv.Accession == CrosslinkItemAccession))
        {
            return "crosslink identification";
        }

        if (string.IsNullOrEmpty(match.PeptideSequence))
        {
            return "peptide reference does not resolve";
        }

        if (match.HasSubstitutionModifications)
        {
            return "peptide carries a substitution modification";
        }

        string sequence = match.PeptideSequence;
        foreach (var modification in match.Modifications)
        {
            if (modification.Location is not int location || location < 0 || location > sequence.Length + 1)
            {
                return $"modification location {modification.Location?.ToString() ?? "(none)"} is not on the peptide";
            }

            var unimod = modification.CvParams.FirstOrDefault(cv => cv.Accession.StartsWith("UNIMOD:", StringComparison.OrdinalIgnoreCase));
            if (unimod == null || !int.TryParse(unimod.Accession["UNIMOD:".Length..], NumberStyles.None, CultureInfo.InvariantCulture, out int unimodId))
            {
                string named = modification.CvParams.Select(cv => cv.Accession).FirstOrDefault(a => a.Length > 0) ?? "(no accession)";
                return $"modification {named} at location {location} has no UNIMOD accession";
            }

            // mzIdentML location: 0 is the N-terminus, 1..length a residue, length + 1 the C-terminus.
            // AllModsOneIsNterminus: 1 is the N-terminus, residue r (zero-based) is r + 2, C-terminus length + 2.
            // Both count the N-terminus as a position, so every key is the location plus one.
            var canonical = location switch
            {
                0 => CanonicalModification.AtNTerminus(unimod.Accession, sequence[0], modification.MonoisotopicMassDelta, unimodId: unimodId),
                _ when location == sequence.Length + 1 => CanonicalModification.AtCTerminus(unimod.Accession, sequence[^1], modification.MonoisotopicMassDelta, unimodId: unimodId),
                _ => CanonicalModification.AtResidue(location - 1, sequence[location - 1], unimod.Accession, modification.MonoisotopicMassDelta, unimodId: unimodId),
            };
            int key = location + 1;

            var resolved = lookup.TryResolve(canonical)?.MzLibModification;
            if (resolved == null)
            {
                return $"{unimod.Accession} at location {location} does not resolve to a Unimod modification";
            }

            if (!mods.TryAdd(key, resolved))
            {
                return $"more than one modification at location {location}";
            }
        }

        return null;
    }

    private static Dictionary<string, double> ScoresOf(MzidSpectrumMatch match)
    {
        var scores = new Dictionary<string, double>();
        var named = match.ItemCvParams.Select(cv => (cv.Name, cv.Value))
            .Concat(match.ItemUserParams.Select(u => (u.Name, u.Value)));
        foreach (var (name, value) in named)
        {
            if (!string.IsNullOrEmpty(name)
                && double.TryParse(value, NumberStyles.Float, CultureInfo.InvariantCulture, out double number))
            {
                scores.TryAdd(name, number);
            }
        }

        return scores;
    }

    private static double? QValueOf(MzidSpectrumMatch match)
    {
        foreach (var cv in match.ItemCvParams)
        {
            if (QValueAccessions.Contains(cv.Accession)
                && double.TryParse(cv.Value, NumberStyles.Float, CultureInfo.InvariantCulture, out double qValue))
            {
                return qValue;
            }
        }

        return null;
    }
}
