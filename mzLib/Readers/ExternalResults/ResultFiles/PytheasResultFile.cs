using Readers.ExternalResults.IndividualResultRecords;
using System.Globalization;
using System.IO;

namespace Readers.ExternalResults.ResultFiles;

/// <summary>
/// A Pytheas match-output file: a '#theoretical_digest'-headed text file whose match lines are
/// grouped under PRECURSOR_ION headers. See <see cref="PytheasResult"/> for the line layout.
/// </summary>
public class PytheasResultFile : ResultFile<PytheasResult>
{
    public override SupportedFileType FileType => SupportedFileType.PytheasResult;
    public override Software Software { get; set; }

    private List<string> _headerLines;

    /// <summary>
    /// The '#' metadata lines at the top of the file, kept so writing round-trips them.
    /// Reading this forces the lazy load, so it is populated immediately after construction.
    /// </summary>
    public List<string> HeaderLines
    {
        get
        {
            EnsureLoaded();
            return _headerLines ??= new List<string>();
        }
    }

    /// <summary>
    /// The number of blank lines preceding each PRECURSOR_ION group, captured during load so that
    /// writing can reproduce the original layout instead of normalizing it to one blank line.
    /// </summary>
    private readonly List<int> _blankLinesBeforeGroups = new();

    #region Header Properties

    private string _theoreticalDigestPath;

    /// <summary>Path of the theoretical digest file, from "#theoretical_digest".</summary>
    public string TheoreticalDigestPath
    {
        get { EnsureLoaded(); return _theoreticalDigestPath; }
        private set => _theoreticalDigestPath = value;
    }

    private string _enzyme;

    /// <summary>Enzyme used for digestion, from "#enzyme".</summary>
    public string Enzyme
    {
        get { EnsureLoaded(); return _enzyme; }
        private set => _enzyme = value;
    }

    private string _msDataPath;

    /// <summary>Path of the MS data file, from "#MS_data".</summary>
    public string MSDataPath
    {
        get { EnsureLoaded(); return _msDataPath; }
        private set => _msDataPath = value;
    }

    private string _isotopicSpecies;

    /// <summary>Isotopic species searched, from "#isotopic_species".</summary>
    public string IsotopicSpecies
    {
        get { EnsureLoaded(); return _isotopicSpecies; }
        private set => _isotopicSpecies = value;
    }

    private double _ms1Ppm;

    /// <summary>MS1 tolerance in ppm, from "#MS1_ppm".</summary>
    public double Ms1Ppm
    {
        get { EnsureLoaded(); return _ms1Ppm; }
        private set => _ms1Ppm = value;
    }

    private double _ms2Ppm;

    /// <summary>MS2 tolerance in ppm, from "#MS2_ppm".</summary>
    public double Ms2Ppm
    {
        get { EnsureLoaded(); return _ms2Ppm; }
        private set => _ms2Ppm = value;
    }

    private double _ms1OffsetPpm;

    /// <summary>MS1 offset tolerance in ppm, from "#MS1_offset_ppm".</summary>
    public double Ms1OffsetPpm
    {
        get { EnsureLoaded(); return _ms1OffsetPpm; }
        private set => _ms1OffsetPpm = value;
    }

    private double _ms2OffsetPpm;

    /// <summary>MS2 offset tolerance in ppm, from "#MS2_offset_ppm".</summary>
    public double Ms2OffsetPpm
    {
        get { EnsureLoaded(); return _ms2OffsetPpm; }
        private set => _ms2OffsetPpm = value;
    }

    private double _ms1MzMinimum;

    /// <summary>MS1 m/z minimum, from "#MS1_mz_minimum".</summary>
    public double Ms1MzMinimum
    {
        get { EnsureLoaded(); return _ms1MzMinimum; }
        private set => _ms1MzMinimum = value;
    }

    private double _ms1MzMaximum;

    /// <summary>MS1 m/z maximum, from "#MS1_mz_maximum".</summary>
    public double Ms1MzMaximum
    {
        get { EnsureLoaded(); return _ms1MzMaximum; }
        private set => _ms1MzMaximum = value;
    }

    private double _ms2MzMinimum;

    /// <summary>MS2 m/z minimum, from "#MS2_mz_minimum".</summary>
    public double Ms2MzMinimum
    {
        get { EnsureLoaded(); return _ms2MzMinimum; }
        private set => _ms2MzMinimum = value;
    }

    private double _ms2MzMaximum;

    /// <summary>MS2 m/z maximum, from "#MS2_mz_maximum".</summary>
    public double Ms2MzMaximum
    {
        get { EnsureLoaded(); return _ms2MzMaximum; }
        private set => _ms2MzMaximum = value;
    }

    private string _ms2AbsPeakIntensity;

    /// <summary>MS2 absolute peak intensity filter, from "#MS2_abs_peak_intensity".</summary>
    public string Ms2AbsPeakIntensity
    {
        get { EnsureLoaded(); return _ms2AbsPeakIntensity; }
        private set => _ms2AbsPeakIntensity = value;
    }

    private string _ms2PeakNumMaximum;

    /// <summary>MS2 peak number maximum, from "#MS2_peak_num_maximum".</summary>
    public string Ms2PeakNumMaximum
    {
        get { EnsureLoaded(); return _ms2PeakNumMaximum; }
        private set => _ms2PeakNumMaximum = value;
    }

    private double _ms2NormintCutoff;

    /// <summary>MS2 normalized intensity cutoff, from "#MS2_normint_cutoff".</summary>
    public double Ms2NormintCutoff
    {
        get { EnsureLoaded(); return _ms2NormintCutoff; }
        private set => _ms2NormintCutoff = value;
    }

    private double _precursorExclusionWindow;

    /// <summary>Precursor exclusion window in Da, from "#precursor_exclusion_window".</summary>
    public double PrecursorExclusionWindow
    {
        get { EnsureLoaded(); return _precursorExclusionWindow; }
        private set => _precursorExclusionWindow = value;
    }

    private double _precursorLossesExclusionWindow;

    /// <summary>Precursor losses exclusion window in Da, from "#precursor_losses_exclusion_window".</summary>
    public double PrecursorLossesExclusionWindow
    {
        get { EnsureLoaded(); return _precursorLossesExclusionWindow; }
        private set => _precursorLossesExclusionWindow = value;
    }

    private double _alpha;

    /// <summary>Pytheas match-scoring alpha parameter, from "#alpha".</summary>
    public double Alpha
    {
        get { EnsureLoaded(); return _alpha; }
        private set => _alpha = value;
    }

    private double _beta;

    /// <summary>Pytheas match-scoring beta parameter, from "#beta".</summary>
    public double Beta
    {
        get { EnsureLoaded(); return _beta; }
        private set => _beta = value;
    }

    private bool _precursorIsotopologues;

    /// <summary>Whether precursor isotopologues were searched, from "#precursor_isotopologues".</summary>
    public bool PrecursorIsotopologues
    {
        get { EnsureLoaded(); return _precursorIsotopologues; }
        private set => _precursorIsotopologues = value;
    }

    private string _matchesHeader;

    /// <summary>Raw column list after "#MATCHES_HEADER:".</summary>
    public string MatchesHeader
    {
        get { EnsureLoaded(); return _matchesHeader; }
        private set => _matchesHeader = value;
    }

    #endregion

    public PytheasResultFile(string filePath) : base(filePath, Software.Pytheas) { }

    /// <summary>Constructor used to initialize from the factory method.</summary>
    public PytheasResultFile() : base() { }

    /// <summary>Forces the lazy file load so callers never observe unpopulated header state.</summary>
    private void EnsureLoaded()
    {
        _ = Results;
    }

    public override void LoadResults()
    {
        _headerLines = new List<string>();
        _blankLinesBeforeGroups.Clear();
        var results = new List<PytheasResult>();
        string precursorIon = null;
        int blankLines = 0;

        foreach (var line in File.ReadLines(FilePath))
        {
            if (line.Length == 0)
            {
                blankLines++;
                continue;
            }

            if (line.StartsWith("#", StringComparison.Ordinal))
            {
                _headerLines.Add(line);
                continue;
            }

            if (line.StartsWith("PRECURSOR_ION=", StringComparison.Ordinal))
            {
                _blankLinesBeforeGroups.Add(blankLines);
                blankLines = 0;
                precursorIon = line.Substring("PRECURSOR_ION=".Length);
                continue;
            }

            results.Add(PytheasResult.Parse(line, precursorIon));
        }

        ParseHeaderLines();
        Results = results;
    }

    /// <summary>
    /// Parses the '#key value' header lines into their typed properties. Lines that don't follow
    /// the shape (e.g. "#MATCHES_HEADER:...", which separates key and value with a colon) are
    /// handled individually; unknown keys are ignored so the reader tolerates future Pytheas versions.
    /// </summary>
    private void ParseHeaderLines()
    {
        foreach (string line in _headerLines)
        {
            string key;
            string value;
            int colon = line.IndexOf(':');
            int space = line.IndexOf(' ');
            if (colon >= 0 && (space < 0 || colon < space))
            {
                key = line.Substring(1, colon - 1);
                value = line.Substring(colon + 1).Trim();
            }
            else if (space >= 0)
            {
                key = line.Substring(1, space - 1);
                value = line.Substring(space + 1).Trim();
            }
            else
            {
                continue;
            }

            switch (key)
            {
                case "theoretical_digest": TheoreticalDigestPath = value; break;
                case "enzyme": Enzyme = value; break;
                case "MS_data": MSDataPath = value; break;
                case "isotopic_species": IsotopicSpecies = value; break;
                case "MS1_ppm": Ms1Ppm = ParseDouble(value); break;
                case "MS2_ppm": Ms2Ppm = ParseDouble(value); break;
                case "MS1_offset_ppm": Ms1OffsetPpm = ParseDouble(value); break;
                case "MS2_offset_ppm": Ms2OffsetPpm = ParseDouble(value); break;
                case "MS1_mz_minimum": Ms1MzMinimum = ParseDouble(value); break;
                case "MS1_mz_maximum": Ms1MzMaximum = ParseDouble(value); break;
                case "MS2_mz_minimum": Ms2MzMinimum = ParseDouble(value); break;
                case "MS2_mz_maximum": Ms2MzMaximum = ParseDouble(value); break;
                case "MS2_abs_peak_intensity": Ms2AbsPeakIntensity = value; break;
                case "MS2_peak_num_maximum": Ms2PeakNumMaximum = value; break;
                case "MS2_normint_cutoff": Ms2NormintCutoff = ParseDouble(value); break;
                case "precursor_exclusion_window": PrecursorExclusionWindow = ParseDouble(value); break;
                case "precursor_losses_exclusion_window": PrecursorLossesExclusionWindow = ParseDouble(value); break;
                case "alpha": Alpha = ParseDouble(value); break;
                case "beta": Beta = ParseDouble(value); break;
                case "precursor_isotopologues": PrecursorIsotopologues = bool.TryParse(value, out bool flag) && flag; break;
                case "MATCHES_HEADER": MatchesHeader = value; break;
            }
        }
    }

    private static double ParseDouble(string value) =>
        double.TryParse(value, NumberStyles.Float, CultureInfo.InvariantCulture, out double result) ? result : 0;

    public override void WriteResults(string outputPath)
    {
        if (!CanRead(outputPath))
            outputPath += FileType.GetFileExtension();

        // Reading Results forces the lazy load, which populates the header lines.
        var results = Results;

        using var writer = new StreamWriter(outputPath);

        foreach (var header in HeaderLines)
            writer.WriteLine(header);

        int groupIndex = 0;
        string lastPrecursorIon = null;
        foreach (var result in results)
        {
            if (result.PrecursorIon != lastPrecursorIon)
            {
                int leadingBlankLines = groupIndex < _blankLinesBeforeGroups.Count
                    ? _blankLinesBeforeGroups[groupIndex]
                    : 1;
                for (int i = 0; i < leadingBlankLines; i++)
                    writer.WriteLine();
                writer.WriteLine("PRECURSOR_ION=" + result.PrecursorIon);
                lastPrecursorIon = result.PrecursorIon;
                groupIndex++;
            }
            writer.WriteLine(result.RawLine);
        }
    }
}