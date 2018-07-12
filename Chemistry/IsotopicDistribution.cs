// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016, 2017 Stefan Solntsev
//
// This file (IsotopicDistribution.cs) is part of Chemistry Library.
//
// Chemistry Library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Chemistry Library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with Chemistry Library. If not, see <http://www.gnu.org/licenses/>.

using System;
using System.Collections.Generic;
using System.Linq;

namespace Chemistry
{
    /// <summary>
    /// Calculates the isotopic distributions of molecules
    /// </summary>
    /// <remarks>
    /// C# version by Derek Bailey 2014
    /// Modified by Stefan Solntsev 2016
    ///
    /// This is a port of software written in C++ and detailed in the following publication:
    ///
    /// Molecular Isotopic Distribution Analysis (MIDAs) with Adjustable Mass Accuracy.
    /// Gelio Alves, Aleksy Y. Ogurtsov, and Yi-Kuo Yu
    /// J. Am. Soc. Mass Spectrom. (2014) 25:57-70
    /// DOI: 10.1007/s13361-013-0733-7
    ///
    /// Please cite that publication if using these algorithms in your own publications.
    ///
    /// Only calculates the fine grained distribution.
    /// </remarks>
    public class IsotopicDistribution
    {
        private const double defaultFineResolution = 0.01;
        private const double defaultMinProbability = 1e-200;
        private const double defaultMolecularWeightResolution = 1e-12;
        private static readonly double[] factorLnArray = new double[50003];
        private static int _factorLnTop = 1;
        private readonly double[] masses;
        private readonly double[] intensities;

        private IsotopicDistribution(int count)
        {
            masses = new double[count];
            intensities = new double[count];
        }

        public IEnumerable<double> Masses { get { return masses; } }
        public IEnumerable<double> Intensities { get { return intensities; } }

        public static IsotopicDistribution GetDistribution(ChemicalFormula formula)
        {
            return GetDistribution(formula, defaultFineResolution, defaultMinProbability, defaultMolecularWeightResolution);
        }

        public static IsotopicDistribution GetDistribution(ChemicalFormula formula, double fineResolution)
        {
            return GetDistribution(formula, fineResolution, defaultMinProbability, defaultMolecularWeightResolution);
        }

        public static IsotopicDistribution GetDistribution(ChemicalFormula formula, double fineResolution, double minProbability)
        {
            return GetDistribution(formula, fineResolution, minProbability, defaultMolecularWeightResolution);
        }

        public static IsotopicDistribution GetDistribution(ChemicalFormula formula, double fineResolution, double minProbability, double molecularWeightResolution)
        {
            var a = GetNewFineAndMergeResolutions(fineResolution);
            var newFineResolution = a.Item1;
            double _mergeFineResolution = a.Item2;
            List<List<Composition>> elementalComposition = new List<List<Composition>>();

            // Get all the unique elements that might have isotopes
            foreach (var elementAndCount in formula.Elements)
            {
                int count = elementAndCount.Value;
                List<Composition> isotopeComposition = new List<Composition>();
                foreach (Isotope isotope in elementAndCount.Key.Isotopes.OrderBy(iso => iso.AtomicMass))
                {
                    Composition c = new Composition
                    {
                        Atoms = count,
                        MolecularWeight = isotope.AtomicMass,
                        Power = isotope.AtomicMass,
                        Probability = isotope.RelativeAbundance
                    };

                    isotopeComposition.Add(c);
                }
                elementalComposition.Add(isotopeComposition);
            }

            foreach (List<Composition> compositions in elementalComposition)
            {
                double sumProb = compositions.Sum(t => t.Probability);
                foreach (Composition composition in compositions)
                {
                    composition.Probability /= sumProb;
                    composition.LogProbability = Math.Log(composition.Probability);
                    composition.Power = Math.Floor(composition.MolecularWeight / molecularWeightResolution + 0.5);
                }
            }
            IsotopicDistribution dist = CalculateFineGrain(elementalComposition, molecularWeightResolution, _mergeFineResolution, newFineResolution, minProbability);

            double additionalMass = 0;
            foreach (var isotopeAndCount in formula.Isotopes)
            {
                additionalMass += isotopeAndCount.Key.AtomicMass * isotopeAndCount.Value;
            }

            for (int i = 0; i < dist.masses.Length; i++)
            {
                dist.masses[i] += additionalMass;
            }

            return dist;
        }

        /// <summary>
        /// Calculates the fineResolution and mergeFineResolution parameters
        /// </summary>
        /// <returns>Tuple of fineResolution and mergeFineResolution</returns>
        private static Tuple<double, double> GetNewFineAndMergeResolutions(double fineResolution)
        {
            return new Tuple<double, double>(fineResolution / 2.0, fineResolution);
        }

        private static List<Polynomial> MergeFinePolynomial(List<Polynomial> tPolynomial, double _mwResolution, double _mergeFineResolution)
        {
            // Sort by mass (i.e. power)
            tPolynomial.Sort((a, b) => a.Power.CompareTo(b.Power));

            int count = tPolynomial.Count;

            for (int k = 1; k <= 9; k++)
            {
                for (int i = 0; i < count; i++)
                {
                    double power = tPolynomial[i].Power;

                    if (double.IsNaN(power))
                        continue;

                    double probability = tPolynomial[i].Probablity;
                    Polynomial tempPolynomial;
                    tempPolynomial.Power = power * probability;
                    tempPolynomial.Probablity = probability;

                    for (int j = i + 1; j < count; j++)
                    {
                        double value = Math.Abs(tPolynomial[i].Power * _mwResolution - tPolynomial[j].Power * _mwResolution);

                        double threshold = (k <= 8) ? k * _mergeFineResolution / 8 : _mergeFineResolution + _mergeFineResolution / 100;

                        // Combine terms if their mass difference (i.e. power difference) is less than some threshold
                        if (value <= threshold)
                        {
                            tempPolynomial.Power = tempPolynomial.Power + tPolynomial[j].Power * tPolynomial[j].Probablity;
                            tempPolynomial.Probablity = tempPolynomial.Probablity + tPolynomial[j].Probablity;
                            tPolynomial[i] = new Polynomial { Power = tempPolynomial.Power / tempPolynomial.Probablity, Probablity = tempPolynomial.Probablity };
                            tPolynomial[j] = new Polynomial { Probablity = double.NaN, Power = double.NaN };
                        }
                        else
                        {
                            break;
                        }
                    }

                    tPolynomial[i] = new Polynomial { Power = tempPolynomial.Power / tempPolynomial.Probablity, Probablity = tempPolynomial.Probablity };
                }
            }

            // return only non-zero terms
            return tPolynomial.Where(poly => !double.IsNaN(poly.Power)).ToList();
        }

        private static List<Polynomial> MultiplyFinePolynomial(List<List<Composition>> elementalComposition, double _fineResolution, double _mwResolution, double _fineMinProb)
        {
            const int nc = 10;
            const int ncAddValue = 1;
            const int nAtoms = 200;
            List<Polynomial> tPolynomial = new List<Polynomial>();

            int n = 0;
            int k;

            foreach (List<Composition> composition in elementalComposition)
                if (composition.Count > 0)
                    n++;

            List<List<Polynomial>> fPolynomial = new List<List<Polynomial>>();
            for (int i = 0; i < n; i++)
                fPolynomial.Add(new List<Polynomial>());

            for (k = 0; k < n; k++)
            {
                tPolynomial.Clear();

                List<Composition> composition = elementalComposition[k];
                int size = composition.Count;
                int atoms = composition[0].Atoms;

                int ncAdd = atoms < nAtoms ? 10 : ncAddValue;

                if (size == 1)
                {
                    double probability = composition[0].Probability;

                    int n1 = (int)(atoms * probability);

                    double prob = FactorLn(atoms) - FactorLn(n1) + n1 * composition[0].LogProbability;
                    prob = Math.Exp(prob);

                    fPolynomial[k].Add(new Polynomial { Power = n1 * composition[0].Power, Probablity = prob });
                }
                else
                {
                    int[] means = new int[size];
                    int[] stds = new int[size];
                    int[] indices = new int[size];

                    double nPolynomialTerms = Math.Log(Math.Pow(2, size));
                    for (int i = 0; i < size; i++)
                    {
                        int n1 = (int)(elementalComposition[k][0].Atoms * elementalComposition[k][i].Probability);
                        int s1 = (int)Math.Ceiling(ncAdd + nc * Math.Sqrt(elementalComposition[k][i].Atoms * elementalComposition[k][i].Probability * (1.0 - elementalComposition[k][i].Probability)));
                        nPolynomialTerms += Math.Log(n1 + s1);

                        means[i] = n1;
                        stds[i] = s1;
                        indices[i] = n1 + s1;
                    }
                    int[] mins = new int[means.Length - 1];
                    int[] maxs = new int[means.Length - 1];
                    indices = new int[means.Length - 1];
                    for (int i = 0; i < means.Length - 1; i++)
                    {
                        var max = Math.Max(0, means[i] - stds[i]);
                        indices[i] = max;
                        mins[i] = max;
                        maxs[i] = means[i] + stds[i];
                    }

                    MultipleFinePolynomialRecursiveHelper(mins, maxs, indices, 0, fPolynomial[k], composition, atoms, _fineMinProb, means[means.Length - 1] + stds[stds.Length - 1]);
                }
            }

            tPolynomial = fPolynomial[0];

            if (k <= 1)
                return tPolynomial;

            List<Polynomial> fgidPolynomial = new List<Polynomial>();
            for (k = 1; k < n; k++)
                MultiplyFineFinalPolynomial(tPolynomial, fPolynomial[k], fgidPolynomial, _fineResolution, _mwResolution, _fineMinProb);

            return tPolynomial;
        }

        private static void MultiplyFineFinalPolynomial(List<Polynomial> tPolynomial, List<Polynomial> fPolynomial, List<Polynomial> fgidPolynomial, double _fineResolution, double _mwResolution, double _fineMinProb)
        {
            int i = tPolynomial.Count;
            int j = fPolynomial.Count;

            if (i == 0 || j == 0)
            {
                return;
            }

            double deltaMass = (_fineResolution / _mwResolution);
            double minProbability = _fineMinProb;

            double minMass = tPolynomial[0].Power + fPolynomial[0].Power;
            double maxMass = tPolynomial[i - 1].Power + fPolynomial[j - 1].Power;

            int maxIndex = (int)(Math.Abs(maxMass - minMass) / deltaMass + 0.5);
            if (maxIndex >= fgidPolynomial.Count)
            {
                j = maxIndex - fgidPolynomial.Count;
                for (i = 0; i <= j; i++)
                {
                    fgidPolynomial.Add(new Polynomial { Probablity = double.NaN, Power = double.NaN });
                }
            }

            for (int t = 0; t < tPolynomial.Count; t++)
            {
                for (int f = 0; f < fPolynomial.Count; f++)
                {
                    double prob = tPolynomial[t].Probablity * fPolynomial[f].Probablity;
                    if (prob <= minProbability)
                    {
                        continue;
                    }

                    double power = tPolynomial[t].Power + fPolynomial[f].Power;
                    var indext = (int)(Math.Abs(power - minMass) / deltaMass + 0.5);

                    Polynomial tempPolynomial = fgidPolynomial[indext];

                    var poww = tempPolynomial.Power;
                    var probb = tempPolynomial.Probablity;
                    if (double.IsNaN(poww) || double.IsNaN(prob))
                    {
                        fgidPolynomial[indext] = new Polynomial { Power = power * prob, Probablity = prob };
                    }
                    else
                    {
                        fgidPolynomial[indext] = new Polynomial { Power = poww + power * prob, Probablity = probb + prob };
                    }
                }
            }

            var index = tPolynomial.Count;
            j = 0;
            for (i = 0; i < fgidPolynomial.Count; i++)
            {
                if (!double.IsNaN(fgidPolynomial[i].Probablity))
                {
                    if (j < index)
                    {
                        tPolynomial[j] = new Polynomial { Power = fgidPolynomial[i].Power / fgidPolynomial[i].Probablity, Probablity = fgidPolynomial[i].Probablity };
                        j++;
                    }
                    else
                        tPolynomial.Add(new Polynomial { Power = fgidPolynomial[i].Power / fgidPolynomial[i].Probablity, Probablity = fgidPolynomial[i].Probablity });
                }

                fgidPolynomial[i] = new Polynomial { Probablity = double.NaN, Power = double.NaN };
            }

            if (j < index)
                tPolynomial.RemoveRange(j, tPolynomial.Count - j);
        }

        private static void MultipleFinePolynomialRecursiveHelper(int[] mins, int[] maxs, int[] indices, int index, IList<Polynomial> fPolynomial, IList<Composition> elementalComposition, int atoms, double minProb, int maxValue)
        {
            for (indices[index] = mins[index]; indices[index] <= maxs[index]; indices[index]++)
            {
                if (index < mins.Length - 1)
                {
                    MultipleFinePolynomialRecursiveHelper(mins, maxs, indices, index + 1, fPolynomial, elementalComposition, atoms, minProb, maxValue);
                }
                else
                {
                    int l = atoms - indices.Sum();
                    if (l < 0 || l > maxValue)
                        continue;

                    double prob = FactorLn(atoms) - FactorLn(l) + l * elementalComposition[elementalComposition.Count - 1].LogProbability;
                    double power = l * elementalComposition[elementalComposition.Count - 1].Power;
                    for (int i = 0; i <= elementalComposition.Count - 2; i++)
                    {
                        int indexValue = indices[i];
                        Composition tComposition = elementalComposition[i];
                        prob -= FactorLn(indexValue);
                        prob += indexValue * tComposition.LogProbability;
                        power += indexValue * tComposition.Power;
                    }

                    prob = Math.Exp(prob);
                    if (prob >= minProb)
                    {
                        Polynomial tPolynomial = new Polynomial { Probablity = prob, Power = power };
                        fPolynomial.Add(tPolynomial);
                    }
                }
            }
        }

        private static double FactorLn(int n)
        {
            if (n <= 1)
            {
                return 0;
            }
            while (_factorLnTop <= n)
            {
                int j = _factorLnTop++;
                factorLnArray[j + 1] = factorLnArray[j] + Math.Log(_factorLnTop);
            }
            return factorLnArray[n];
        }

        private static IsotopicDistribution CalculateFineGrain(List<List<Composition>> elementalComposition, double _mwResolution, double _mergeFineResolution, double _fineResolution, double _fineMinProb)
        {
            List<Polynomial> fPolynomial = MultiplyFinePolynomial(elementalComposition, _fineResolution, _mwResolution, _fineMinProb);
            fPolynomial = MergeFinePolynomial(fPolynomial, _mwResolution, _mergeFineResolution);

            // Convert polynomial to spectrum
            int count = fPolynomial.Count;
            IsotopicDistribution dist = new IsotopicDistribution(count);
            double totalProbability = 0;
            double basePeak = 0;
            int i = 0;
            foreach (Polynomial polynomial in fPolynomial)
            {
                totalProbability += polynomial.Probablity;
                if (polynomial.Probablity > basePeak)
                {
                    basePeak = polynomial.Probablity;
                }
                dist.masses[i] = polynomial.Power * _mwResolution;
                dist.intensities[i] = polynomial.Probablity;
                i++;
            }
            return dist;
        }

        private struct Polynomial
        {
            public double Power;
            public double Probablity;
        }

        private class Composition
        {
            public double Power;
            public double Probability;
            public double LogProbability;
            public double MolecularWeight;
            public int Atoms;
        }
    }
}