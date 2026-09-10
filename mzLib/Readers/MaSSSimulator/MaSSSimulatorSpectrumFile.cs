using System.Globalization;

namespace Readers.MaSSSimulator;

/// <summary>
/// Reads MaSS-Simulator's H/S/Z/peak text output and converts it to MGF.
/// MaSS-Simulator writes truth to a separate <c>peptides.rst</c> file, so
/// provenance is loaded separately and never inferred from peak data.
/// </summary>
public sealed class MaSSSimulatorSpectrumFile : ResultFile<MaSSSimulatorSpectrum>
{
    public override SupportedFileType FileType => SupportedFileType.MaSSSimulatorSpectra;
    public override Software Software { get; set; } = Software.Unspecified;

    public MaSSSimulatorSpectrumFile() : base() { }
    public MaSSSimulatorSpectrumFile(string filePath) : base(filePath) { }
    public MaSSSimulatorSpectrumFile(IEnumerable<MaSSSimulatorSpectrum> spectra) : base() => Results = spectra.ToList();

    public override void LoadResults()
    {
        var spectra = new List<MaSSSimulatorSpectrum>();
        MaSSSimulatorSpectrum? current = null;

        foreach (var rawLine in File.ReadLines(FilePath))
        {
            var line = rawLine.Trim();
            if (line.Length == 0 || line.StartsWith("H\t", StringComparison.Ordinal))
                continue;

            var fields = line.Split('\t');
            switch (fields[0])
            {
                case "S":
                    if (fields.Length < 4)
                        throw new FormatException($"Invalid MaSS-Simulator spectrum line: {rawLine}");
                    current = new MaSSSimulatorSpectrum
                    {
                        ScanNumber = ParseInt(fields[1], rawLine),
                        PrecursorMz = ParseDouble(fields[3], rawLine)
                    };
                    spectra.Add(current);
                    break;
                case "Z":
                    if (current is null || fields.Length < 3)
                        throw new FormatException($"MaSS-Simulator Z line without an S line: {rawLine}");
                    current.Charge = ParseInt(fields[1], rawLine);
                    current.PrecursorMass = ParseDouble(fields[2], rawLine);
                    break;
                default:
                    if (current is null)
                        throw new FormatException($"Peak line without an S line: {rawLine}");
                    var peak = line.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
                    if (peak.Length != 2)
                        throw new FormatException($"Invalid MaSS-Simulator peak line: {rawLine}");
                    current.Peaks.Add((ParseDouble(peak[0], rawLine), ParseDouble(peak[1], rawLine)));
                    break;
            }
        }

        Results = spectra;
    }

    public override void WriteResults(string outputPath)
    {
        using var writer = new StreamWriter(outputPath);
        writer.WriteLine("H\tCreationDate\t");
        foreach (var spectrum in Results)
        {
            writer.WriteLine($"S\t{spectrum.ScanNumber}\t{spectrum.ScanNumber}\t{Format(spectrum.PrecursorMz)}");
            if (spectrum.PrecursorMass.HasValue)
                writer.WriteLine($"Z\t{spectrum.Charge}\t{Format(spectrum.PrecursorMass.Value)}");
            foreach (var peak in spectrum.Peaks)
                writer.WriteLine($"{Format(peak.Mz)} {Format(peak.Intensity)}");
        }
    }

    /// <summary>Writes spectra as MS2-only Mascot Generic Format.</summary>
    public void WriteMgf(string outputPath)
    {
        using var writer = new StreamWriter(outputPath);
        foreach (var spectrum in Results.Where(s => s.Peaks.Count > 0))
        {
            writer.WriteLine("BEGIN IONS");
            writer.WriteLine($"TITLE=MaSS-Simulator scan {spectrum.ScanNumber}");
            writer.WriteLine($"PEPMASS={Format(spectrum.PrecursorMz)}");
            if (spectrum.Charge != 0)
                writer.WriteLine($"CHARGE={Math.Abs(spectrum.Charge)}{(spectrum.Charge < 0 ? '-' : '+')}");
            writer.WriteLine($"SCANS={spectrum.ScanNumber}");
            foreach (var peak in spectrum.Peaks)
                writer.WriteLine($"{Format(peak.Mz)} {Format(peak.Intensity)}");
            writer.WriteLine("END IONS");
        }
    }

    public void ApplyTruth(string rstPath)
    {
        foreach (var line in File.ReadLines(rstPath))
        {
            var separator = line.IndexOf(" peptide:", StringComparison.Ordinal);
            if (!line.StartsWith("scan:", StringComparison.Ordinal) || separator < 0)
                continue;
            var scan = int.Parse(line[5..separator], CultureInfo.InvariantCulture);
            var peptide = line[(separator + 9)..];
            var spectrum = Results.FirstOrDefault(s => s.ScanNumber == scan);
            if (spectrum is not null)
                spectrum.PeptideSequence = peptide;
        }
    }

    private static int ParseInt(string value, string line) =>
        int.TryParse(value, NumberStyles.Integer, CultureInfo.InvariantCulture, out var parsed)
            ? parsed : throw new FormatException($"Invalid integer in MaSS-Simulator line: {line}");

    private static double ParseDouble(string value, string line) =>
        double.TryParse(value, NumberStyles.Float, CultureInfo.InvariantCulture, out var parsed)
            ? parsed : throw new FormatException($"Invalid number in MaSS-Simulator line: {line}");

    private static string Format(double value) => value.ToString("0.######", CultureInfo.InvariantCulture);
}
