namespace MassSpectrometry;

/// <summary>
/// Represents the average residue and is used for most abundant isotopic peak to monoisotopic peak difference
/// </summary>
/// <remarks>CAREFUL: This is called a lot during <see cref="ClassicDeconvolutionAlgorithm"/> and care should be taken to precalculate as many values as possible when implementing a new type</remarks>
public abstract class AverageResidue
{
    protected const double FineRes = 0.125;
    protected const double MinRes = 1e-8;
    protected static readonly int NumAveraginesToGenerate = 1500;
    public abstract int GetMostIntenseMassIndex(double testMass);
    public abstract double[] GetAllTheoreticalMasses(int index);
    public abstract double[] GetAllTheoreticalIntensities(int index);
    public abstract double GetDiffToMonoisotopic(int index);
}