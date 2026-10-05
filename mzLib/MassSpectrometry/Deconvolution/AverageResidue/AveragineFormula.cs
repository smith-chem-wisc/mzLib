#nullable enable
using Chemistry;
using System;
using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry;

/// <summary>
/// Stands in averagine atoms for the part of a species' mass that has no known chemical formula, so that a
/// theoretical isotopic envelope can still be computed for it.
/// </summary>
/// <remarks>
/// The per-residue compositions come from <see cref="Averagine.GetAverageChemicalFormula"/> (amino acids) and
/// <see cref="OxyriboAveragine.GetAverageChemicalFormula"/> (ribonucleotides). FlashLFQ uses this for identifications
/// whose modifications carry a mass but no formula.
/// </remarks>
public static class AveragineFormula
{
    /// <summary>
    /// Adds <paramref name="mass"/> daltons of averagine to <paramref name="formula"/>: the per-residue
    /// <paramref name="averageComposition"/> is scaled by the number of averagine residues in <paramref name="mass"/>,
    /// measured with element average masses, and each element is rounded to the nearest whole atom.
    /// </summary>
    /// <remarks>
    /// A negative <paramref name="mass"/> subtracts atoms, and element counts are not clamped at zero.
    /// </remarks>
    public static void AddAveragine(ChemicalFormula formula, double mass, IReadOnlyDictionary<char, double> averageComposition)
    {
        ArgumentNullException.ThrowIfNull(formula);
        ArgumentNullException.ThrowIfNull(averageComposition);

        double averagineMass = averageComposition
            .Sum(kvp => PeriodicTable.GetElement(kvp.Key.ToString()).AverageMass * kvp.Value);
        double averagines = mass / averagineMass;
        foreach (var (element, countPerAveragine) in averageComposition)
        {
            formula.Add(element.ToString(), (int)Math.Round(averagines * countPerAveragine, 0));
        }
    }

    /// <summary>
    /// Returns the theoretical isotopic envelope of a species whose exact monoisotopic mass is
    /// <paramref name="monoisotopicMass"/> and whose known <paramref name="formula"/> may not account for all of it, such
    /// as a peptide carrying a modification that has a mass but no formula. When the two masses differ by more than
    /// <paramref name="averagineThreshold"/> daltons in either direction, the difference is filled with averagine
    /// (<see cref="AddAveragine"/>) on a copy of <paramref name="formula"/>. Every mass is then shifted so the
    /// monoisotopic peak sits exactly at <paramref name="monoisotopicMass"/>.
    /// </summary>
    /// <param name="formula">The known part of the species. It is not modified.</param>
    /// <param name="monoisotopicMass">The species' exact monoisotopic mass.</param>
    /// <param name="averageComposition">The averagine composition to fill with.</param>
    /// <param name="fineResolution">Passed to <see cref="IsotopicDistribution.GetDistribution(ChemicalFormula, double, double)"/>.</param>
    /// <param name="minProbability">Passed to <see cref="IsotopicDistribution.GetDistribution(ChemicalFormula, double, double)"/>.</param>
    /// <param name="averagineThreshold">The largest mass difference, in daltons, left unfilled. Pass
    /// <see cref="double.PositiveInfinity"/> to never fill.</param>
    /// <returns>The masses in ascending order and their intensities, as
    /// <see cref="IsotopicDistribution.GetDistribution(ChemicalFormula, double, double)"/> returns them.</returns>
    public static (double[] Masses, double[] Intensities) GetAnchoredDistribution(ChemicalFormula formula, double monoisotopicMass,
        IReadOnlyDictionary<char, double> averageComposition, double fineResolution, double minProbability, double averagineThreshold = 20)
    {
        ArgumentNullException.ThrowIfNull(formula);
        ArgumentNullException.ThrowIfNull(averageComposition);

        var filled = new ChemicalFormula(formula);
        double massDiff = monoisotopicMass - filled.MonoisotopicMass;
        if (Math.Abs(massDiff) > averagineThreshold)
        {
            AddAveragine(filled, massDiff, averageComposition);
        }

        var distribution = IsotopicDistribution.GetDistribution(filled, fineResolution, minProbability);
        double[] masses = distribution.Masses.ToArray();
        double[] intensities = distribution.Intensities.ToArray();
        for (int i = 0; i < masses.Length; i++)
        {
            masses[i] += (monoisotopicMass - filled.MonoisotopicMass);
        }

        return (masses, intensities);
    }
}
