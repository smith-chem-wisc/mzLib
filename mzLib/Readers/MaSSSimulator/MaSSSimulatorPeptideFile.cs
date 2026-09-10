using Omics;
using Omics.Digestion;
using Omics.Modifications;

namespace Readers.MaSSSimulator;

/// <summary>
/// Reads and writes the peptide-list format consumed by MaSS-Simulator's SimSpec.
/// The simulator requires a header line followed by one peptide sequence per line.
/// Provenance is retained in memory and can be written as a separate tab-delimited
/// truth table; it is deliberately not placed in the simulator input file.
/// </summary>
public sealed class MaSSSimulatorPeptideFile : ResultFile<MaSSSimulatorPeptide>
{
    public override SupportedFileType FileType => SupportedFileType.MaSSSimulatorPeptides;
    public override Software Software { get; set; } = Software.Unspecified;

    public MaSSSimulatorPeptideFile(string filePath) : base(filePath) { }
    public MaSSSimulatorPeptideFile(IEnumerable<MaSSSimulatorPeptide> peptides) : base() => Results = peptides.ToList();

    public override void LoadResults()
    {
        Results = File.ReadLines(FilePath)
            .Skip(1)
            .Where(line => !string.IsNullOrWhiteSpace(line))
            .Select(line => new MaSSSimulatorPeptide(line.Trim()))
            .ToList();
    }

    public override void WriteResults(string outputPath)
    {
        using var writer = new StreamWriter(outputPath);
        writer.WriteLine("peptide");
        foreach (var peptide in Results)
            writer.WriteLine(peptide.Sequence);
    }

    /// <summary>Creates simulator input from mzLib peptide objects.</summary>
    public static MaSSSimulatorPeptideFile FromPeptides(IEnumerable<IBioPolymerWithSetMods> peptides) =>
        new(peptides.Select(p => new MaSSSimulatorPeptide(ToSimulatorSequence(p), p.Parent?.Accession)));

    /// <summary>
    /// Digests mzLib proteins and creates simulator input. Fixed and variable
    /// modifications are passed through mzLib's normal digestion API.
    /// </summary>
    public static MaSSSimulatorPeptideFile FromProteins(
        IEnumerable<IBioPolymer> proteins,
        IDigestionParams digestionParams,
        List<Modification>? fixedModifications = null,
        List<Modification>? variableModifications = null) =>
        FromPeptides(proteins.SelectMany(p => p.Digest(
            digestionParams,
            fixedModifications ?? [],
            variableModifications ?? [])));

    /// <summary>Writes the source accession alongside each simulator peptide.</summary>
    public void WriteTruthTable(string outputPath)
    {
        using var writer = new StreamWriter(outputPath);
        writer.WriteLine("scan\tpeptide\taccession");
        for (var i = 0; i < Results.Count; i++)
            writer.WriteLine($"{i + 1}\t{Results[i].Sequence}\t{Results[i].Accession ?? ""}");
    }

    private static string ToSimulatorSequence(IBioPolymerWithSetMods peptide)
    {
        if (peptide.AllModsOneIsNterminus.Count > 0)
            throw new NotSupportedException("MaSS-Simulator adapter currently supports unmodified peptides only; PTM mass conversion is planned.");
        return peptide.BaseSequence;
    }
}
