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

    /// <summary>
    /// The number of daltons to add to a monoisotopic mass to reach the most abundant isotopic peak,
    /// according to this average-residue model. Use this to convert an exact candidate monoisotopic
    /// mass into the theoretical most-abundant mass for matching against an observed most-abundant
    /// peak (Strategy B). The model only supplies the (mass-dependent, roughly integer-neutron)
    /// offset; the candidate monoisotopic mass stays exact.
    /// </summary>
    /// <param name="monoisotopicMass">A monoisotopic mass, in daltons.</param>
    public double GetMostAbundantOffset(double monoisotopicMass)
        => GetDiffToMonoisotopic(GetMostIntenseMassIndex(monoisotopicMass));

    /// <summary>
    /// The number of daltons to add to a monoisotopic mass to reach the intensity-weighted average
    /// (centroid) mass of the isotopic envelope, according to this average-residue model. Used for
    /// isotopically unresolved (high-mass) species, where the observed precursor is a centroid rather
    /// than a resolved most-abundant peak.
    /// </summary>
    /// <param name="monoisotopicMass">A monoisotopic mass, in daltons.</param>
    public double GetAverageOffset(double monoisotopicMass)
    {
        int index = GetMostIntenseMassIndex(monoisotopicMass);
        double[] masses = GetAllTheoreticalMasses(index);
        double[] intensities = GetAllTheoreticalIntensities(index);

        double totalIntensity = 0;
        double weightedMassSum = 0;
        for (int i = 0; i < masses.Length; i++)
        {
            totalIntensity += intensities[i];
            weightedMassSum += masses[i] * intensities[i];
        }
        double averageMass = totalIntensity > 0 ? weightedMassSum / totalIntensity : masses[0];

        // masses/intensities are sorted by intensity (most intense first); the model's monoisotopic
        // mass is the most-intense mass minus its diff-to-monoisotopic.
        double monoisotopicModelMass = masses[0] - GetDiffToMonoisotopic(index);
        return averageMass - monoisotopicModelMass;
    }
}