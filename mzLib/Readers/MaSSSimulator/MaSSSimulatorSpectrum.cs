namespace Readers.MaSSSimulator;

/// <summary>One MS/MS spectrum emitted by MaSS-Simulator.</summary>
public sealed class MaSSSimulatorSpectrum
{
    public int ScanNumber { get; set; }
    public double PrecursorMz { get; set; }
    public int Charge { get; set; }
    public double? PrecursorMass { get; set; }
    public List<(double Mz, double Intensity)> Peaks { get; set; } = [];
    public string? PeptideSequence { get; set; }
    public string? Accession { get; set; }
}
