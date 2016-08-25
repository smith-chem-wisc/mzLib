// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
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
using System.Collections.ObjectModel;
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

        private double[] masses;
        private double[] intensities;

        public ReadOnlyCollection<double> Masses
        {
            get { return new ReadOnlyCollection<double>(masses); }
        }

        public ReadOnlyCollection<double> Intensities
        {
            get { return new ReadOnlyCollection<double>(intensities); }
        }

        public IsotopicDistribution(ChemicalFormula formula) : this(ValidateFormulaForIsotopologueComputation(formula), defaultFineResolution, defaultMinProbability, defaultMolecularWeightResolution)
        {
        }

        public IsotopicDistribution(ChemicalFormula formula, double fineResolution) : this(ValidateFormulaForIsotopologueComputation(formula), fineResolution, defaultMinProbability, defaultMolecularWeightResolution)
        {
        }

        public IsotopicDistribution(ChemicalFormula formula, double fineResolution, double minProbability) : this(ValidateFormulaForIsotopologueComputation(formula), fineResolution, minProbability, defaultMolecularWeightResolution)
        {
        }

        public IsotopicDistribution(ChemicalFormula formula, double fineResolution, double minProbability, double molecularWeightResolution)
        {
            ValidateFormulaForIsotopologueComputation(formula);
            double monoisotopicMass = formula.MonoisotopicMass;
            var a = GetNewFineAndMergeResolutions(monoisotopicMass, fineResolution);
            fineResolution = a.Item1;
            double _mergeFineResolution = a.Item2;
            List<List<Composition>> elementalComposition = new List<List<Composition>>();

            // Get all the unique elements that might have isotopes
            foreach (var elementAndCount in formula.elements)
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
            CalculateFineGrain(elementalComposition, molecularWeightResolution, _mergeFineResolution, fineResolution, minProbability);

            double additionalMass = 0;
            foreach (var isotopeAndCount in formula.isotopes)
                additionalMass += isotopeAndCount.Key.AtomicMass * isotopeAndCount.Value;

            for (int i = 0; i < Masses.Count(); i++)
                masses[i] += additionalMass;
        }

        private static ChemicalFormula ValidateFormulaForIsotopologueComputation(ChemicalFormula formula)
        {
            if (formula == null)
                throw new ArgumentNullException("formula", "Cannot compute isotopic distribution for a null formula");
            return formula;
        }

        /// <summary>
        /// Calculates the fineResolution and mergeFineResolution parameters
        /// </summary>
        /// <remarks>
        /// Takes your guess for fine resolution, and messes it up
        /// If smaller than 1e-4, set to 1e-4
        /// If smaller than 1e-3, set to 1e-3
        /// If smaller than 1e-2, set to 1e-2
        /// If between 1e-2 and 1, set to 1e-2, but merge fine resolution set to original fine resolution
        /// If bigger than 1, set to 0.9
        /// DIVIDE FINE RESOLUTION BY 2 AT THE END
        /// </remarks>
        /// <returns>Tuple of fineResolution and mergeFineResolution</returns>
        private static Tuple<double, double> GetNewFineAndMergeResolutions(double monoisotopicMass, double fineResolution)
        {
            return new Tuple<double, double>(fineResolution / 2.0, fineResolution);
        }

        private void CalculateFineGrain(List<List<Composition>> elementalComposition, double _mwResolution, double _mergeFineResolution, double _fineResolution, double _fineMinProb)
        {
            List<Polynomial> fPolynomial = MultiplyFinePolynomial(elementalComposition, _fineResolution, _mwResolution, _fineMinProb);
            fPolynomial = MergeFinePolynomial(fPolynomial, _mwResolution, _mergeFineResolution);

            // Convert polynomial to spectrum
            int count = fPolynomial.Count;
            masses = new double[count];
            intensities = new double[count];
            double totalProbability = 0;
            double basePeak = 0;
            int i = 0;
            foreach (Polynomial polynomial in fPolynomial)
            {
                totalProbability += polynomial.Probablity;
                if (polynomial.Probablity > basePeak)
                    basePeak = polynomial.Probablity;
                masses[i] = polynomial.Power * _mwResolution;
                intensities[i] = polynomial.Probablity;
                i++;
            }
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

                    if (power == 0)
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
                            tPolynomial[j] = new Polynomial();
                        }
                        else
                            break;
                    }

                    tPolynomial[i] = new Polynomial { Power = tempPolynomial.Power / tempPolynomial.Probablity, Probablity = tempPolynomial.Probablity };
                }
            }

            // return only non-zero terms
            return tPolynomial.Where(poly => poly.Power != 0).ToList();
        }

        private static List<Polynomial> MultiplyFinePolynomial(List<List<Composition>> elementalComposition, double _fineResolution, double _mwResolution, double _fineMinProb)
        {
            const int nc = 10;
            const int ncAddValue = 1;
            const int nAtoms = 200;
            List<Polynomial> tPolynomial = new List<Polynomial>();

            int n = 0;
            int k = 0;

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
                        indices[i] = mins[i] = Math.Max(0, means[i] - stds[i]);
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
                return;

            double deltaMass = (_fineResolution / _mwResolution);
            double minProbability = _fineMinProb;

            int index = 0;
            double minMass = tPolynomial[0].Power + fPolynomial[0].Power;
            double maxMass = tPolynomial[i - 1].Power + fPolynomial[j - 1].Power;

            int maxIndex = (int)(Math.Abs(maxMass - minMass) / deltaMass + 0.5);
            if (maxIndex >= fgidPolynomial.Count)
            {
                j = maxIndex - fgidPolynomial.Count;
                for (i = 0; i <= j; i++)
                    fgidPolynomial.Add(new Polynomial());
            }

            for (int t = 0; t < tPolynomial.Count; t++)
            {
                for (int f = 0; f < fPolynomial.Count; f++)
                {
                    double prob = tPolynomial[t].Probablity * fPolynomial[f].Probablity;
                    if (prob <= minProbability)
                        continue;

                    double power = tPolynomial[t].Power + fPolynomial[f].Power;
                    index = (int)(Math.Abs(power - minMass) / deltaMass + 0.5);

                    Polynomial tempPolynomial = fgidPolynomial[index];

                    fgidPolynomial[index] = new Polynomial { Power = tempPolynomial.Power + power * prob, Probablity = tempPolynomial.Probablity + prob };
                }
            }

            index = tPolynomial.Count;
            j = 0;
            for (i = 0; i < fgidPolynomial.Count; i++)
            {
                if (fgidPolynomial[i].Probablity != 0)
                {
                    if (j < index)
                    {
                        tPolynomial[j] = new Polynomial { Power = fgidPolynomial[i].Power / fgidPolynomial[i].Probablity, Probablity = fgidPolynomial[i].Probablity };
                        j++;
                    }
                    else
                        tPolynomial.Add(new Polynomial { Power = fgidPolynomial[i].Power / fgidPolynomial[i].Probablity, Probablity = fgidPolynomial[i].Probablity });
                }

                fgidPolynomial[i] = new Polynomial();
            }

            if (j < index)
                tPolynomial.RemoveRange(j, tPolynomial.Count - j);
        }

        private static void MultipleFinePolynomialRecursiveHelper(int[] mins, int[] maxs, int[] indices, int index, IList<Polynomial> fPolynomial, IList<Composition> elementalComposition, int atoms, double minProb, int maxValue)
        {
            for (indices[index] = mins[index]; indices[index] <= maxs[index]; indices[index]++)
            {
                if (index < mins.Length - 1)
                    MultipleFinePolynomialRecursiveHelper(mins, maxs, indices, index + 1, fPolynomial, elementalComposition, atoms, minProb, maxValue);
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

        private static readonly double[] factorLnArray = new double[50003];

        private static int _factorLnTop = 1;

        private static double FactorLn(int n)
        {
            if (n <= 1)
                return 0;

            //if (n > 50000)
            //    return n * Math.Log(n) - n + 0.5 * Math.Log(6.28318530717959 * n) + 0.08333333333333 / n - 0.00333333333333 / (n * n * n);

            while (_factorLnTop <= n)
            {
                int j = _factorLnTop++;
                factorLnArray[j + 1] = factorLnArray[j] + Math.Log(_factorLnTop);
            }
            return factorLnArray[n];
        }

        private class Composition
        {
            public double Power;
            public double Probability;
            public double LogProbability;
            public double MolecularWeight;
            public int Atoms;
        }

        // TODO: Benchmark class vs struct here...
        private struct Polynomial
        {
            public double Power;
            public double Probablity;
        }
    }
}