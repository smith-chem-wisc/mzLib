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

    /// <summary>The '#' metadata lines at the top of the file, kept so writing round-trips them.</summary>
    public List<string> HeaderLines { get; } = new();

    #region Header Properties

    /// <summary>Path of the theoretical digest file, from "#theoretical_digest".</summary>
    public string TheoreticalDigestPath { get; private set; }

    /// <summary>Enzyme used for digestion, from "#enzyme".</summary>
    public string Enzyme { get; private set; }

    /// <summary>Path of the MS data file, from "#MS_data".</summary>
    public string MSDataPath { get; private set; }

    /// <summary>Isotopic species searched, from "#isotopic_species".</summary>
    public string IsotopicSpecies { get; private set; }

    /// <summary>MS1 tolerance in ppm, from "#MS1_ppm".</summary>
    public double Ms1Ppm { get; private set; }

    /// <summary>MS2 tolerance in ppm, from "#MS2_ppm".</summary>
    public double Ms2Ppm { get; private set; }

    /// <summary>MS1 offset tolerance in ppm, from "#MS1_offset_ppm".</summary>
    public double Ms1OffsetPpm { get; private set; }

    /// <summary>MS2 offset tolerance in ppm, from "#MS2_offset_ppm".</summary>
    public double Ms2OffsetPpm { get; private set; }

    /// <summary>MS1 m/z minimum, from "#MS1_mz_minimum".</summary>
    public double Ms1MzMinimum { get; private set; }

    /// <summary>MS1 m/z maximum, from "#MS1_mz_maximum".</summary>
    public double Ms1MzMaximum { get; private set; }

    /// <summary>MS2 m/z minimum, from "#MS2_mz_minimum".</summary>
    public double Ms2MzMinimum { get; private set; }

    /// <summary>MS2 m/z maximum, from "#MS2_mz_maximum".</summary>
    public double Ms2MzMaximum { get; private set; }

    /// <summary>MS2 absolute peak intensity filter, from "#MS2_abs_peak_intensity".</summary>
    public string Ms2AbsPeakIntensity { get; private set; }

    /// <summary>MS2 peak number maximum, from "#MS2_peak_num_maximum".</summary>
    public string Ms2PeakNumMaximum { get; private set; }

    /// <summary>MS2 normalized intensity cutoff, from "#MS2_normint_cutoff".</summary>
    public double Ms2NormintCutoff { get; private set; }

    /// <summary>Precursor exclusion window in Da, from "#precursor_exclusion_window".</summary>
    public double PrecursorExclusionWindow { get; private set; }

    /// <summary>Precursor losses exclusion window in Da, from "#precursor_losses_exclusion_window".</summary>
    public double PrecursorLossesExclusionWindow { get; private set; }

    /// <summary>Pytheas match-scoring alpha parameter, from "#alpha".</summary>
    public double Alpha { get; private set; }

    /// <summary>Pytheas match-scoring beta parameter, from "#beta".</summary>
    public double Beta { get; private set; }

    /// <summary>Whether precursor isotopologues were searched, from "#precursor_isotopologues".</summary>
    public bool PrecursorIsotopologues { get; private set; }

    /// <summary>Raw column list after "#MATCHES_HEADER:".</summary>
    public string MatchesHeader { get; private set; }

    #endregion

    public PytheasResultFile(string filePath) : base(filePath, Software.Pytheas) { }

    /// <summary>Constructor used to initialize from the factory method.</summary>
    public PytheasResultFile() : base() { }

    public override void LoadResults()
    {
        var results = new List<PytheasResult>();
        string precursorIon = null;

        foreach (var line in File.ReadLines(FilePath))
        {
            if (line.Length == 0)
                continue;

            if (line.StartsWith("#", StringComparison.Ordinal))
            {
                HeaderLines.Add(line);
                continue;
            }

            if (line.StartsWith("PRECURSOR_ION=", StringComparison.Ordinal))
            {
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
        foreach (string line in HeaderLines)
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

        // Reading Results forces the lazy load, which populates HeaderLines.
        var results = Results;

        using var writer = new StreamWriter(outputPath);

        foreach (var header in HeaderLines)
            writer.WriteLine(header);

        string lastPrecursorIon = null;
        foreach (var result in results)
        {
            if (result.PrecursorIon != lastPrecursorIon)
            {
                writer.WriteLine();
                writer.WriteLine("PRECURSOR_ION=" + result.PrecursorIon);
                lastPrecursorIon = result.PrecursorIon;
            }
            writer.WriteLine(result.RawLine);
        }
    }
}