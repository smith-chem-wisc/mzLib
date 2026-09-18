namespace Readers.MaSSSimulator;

/// <summary>A peptide line and optional mzLib provenance for MaSS-Simulator input.</summary>
public sealed record MaSSSimulatorPeptide(string Sequence, string? Accession = null);
