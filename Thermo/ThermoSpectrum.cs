// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work Copyright 2016 Stefan Solntsev
//
// This file (ThermoSpectrum.cs) is part of MassSpecFiles.
//
// MassSpecFiles is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MassSpecFiles is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with MassSpecFiles. If not, see <http://www.gnu.org/licenses/>.

using Chemistry;
using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;

namespace IO.Thermo
{
    /// <summary>
    /// A high resolution spectra from a Thermo raw file
    /// </summary>
    [Serializable]
    public sealed class ThermoSpectrum : MzSpectrum<ThermoMzPeak>
    {

        #region Private Fields

        private readonly double[] _noises;

        private readonly double[] _resolutions;

        private readonly int[] _charges;

        #endregion Private Fields

        #region Public Constructors

        public ThermoSpectrum(double[] mz, double[] intensity, double[] noise, int[] charge, double[] resolutions, bool shouldCopy)
            : base(mz, intensity, shouldCopy)
        {
            if (!shouldCopy)
            {
                _noises = noise;
                _resolutions = resolutions;
                _charges = charge;
            }
            else
            {
                if (noise != null)
                {
                    _noises = new double[noise.Length];
                    Array.Copy(noise, _noises, noise.Length);
                }
                if (resolutions != null)
                {
                    _resolutions = new double[resolutions.Length];
                    Array.Copy(resolutions, _resolutions, resolutions.Length);
                }
                if (charge != null)
                {
                    _charges = new int[charge.Length];
                    Array.Copy(charge, _charges, charge.Length);
                }
            }
        }

        public ThermoSpectrum(ThermoSpectrum thermoSpectrum)
            : this(thermoSpectrum.XArray, thermoSpectrum.YArray, thermoSpectrum._noises, thermoSpectrum._charges, thermoSpectrum._resolutions, true)
        {
        }

        #endregion Public Constructors

        #region Internal Constructors

        internal ThermoSpectrum(double[,] peakData)
                            : base(peakData)
        {
            int arrayLength = peakData.GetLength(1);
            int depthLength = peakData.GetLength(0);
            if (depthLength <= 2)
                return;
            _noises = new double[Size];
            _resolutions = new double[Size];
            var charges = new double[Size];

            Buffer.BlockCopy(peakData, sizeof(double) * arrayLength * (int)RawLabelDataColumn.Resolution, _resolutions, 0, sizeof(double) * Size);
            Buffer.BlockCopy(peakData, sizeof(double) * arrayLength * (int)RawLabelDataColumn.NoiseLevel, _noises, 0, sizeof(double) * Size);
            Buffer.BlockCopy(peakData, sizeof(double) * arrayLength * (int)RawLabelDataColumn.Charge, charges, 0, sizeof(double) * Size);

            _charges = new int[Size];
            for (int i = 0; i < Size; i++)
            {
                _charges[i] = (int)charges[i];
            }
        }

        #endregion Internal Constructors

        #region Private Enums

        private enum RawLabelDataColumn
        {
            MZ = 0,
            Intensity = 1,
            Resolution = 2,
            NoiseBaseline = 3,
            NoiseLevel = 4,
            Charge = 5
        }

        #endregion Private Enums

        #region Public Methods

        public double GetSignalToNoise(int index)
        {
            if (_noises == null)
                return double.NaN;
            double noise = _noises[index];
            return YArray[index] / noise;
        }

        public double[] GetNoises()
        {
            return _noises;
        }

        public double[] GetResolutions()
        {
            return _resolutions;
        }

        public int[] GetCharges()
        {
            return _charges;
        }

        public override double[,] CopyTo2DArray()
        {
            double[,] data = new double[5, Size];
            const int size = sizeof(double);
            int bytesToCopy = size * Size;
            Buffer.BlockCopy(XArray, 0, data, 0, bytesToCopy);
            Buffer.BlockCopy(YArray, 0, data, bytesToCopy, bytesToCopy);
            Buffer.BlockCopy(_resolutions, 0, data, 2 * bytesToCopy, bytesToCopy);
            Buffer.BlockCopy(_noises, 0, data, 3 * bytesToCopy, bytesToCopy);

            double[] charges = new double[Size];
            for (int i = 0; i < Size; i++)
                charges[i] = _charges[i];

            Buffer.BlockCopy(charges, 0, data, 4 * bytesToCopy, bytesToCopy);
            return data;
        }

        #endregion Public Methods

        #region Protected Methods

        protected override ThermoMzPeak GeneratePeak(int index)
        {
            return peakList[index] = _charges == null ?
                new ThermoMzPeak(XArray[index], YArray[index]) :
                new ThermoMzPeak(XArray[index], YArray[index], _charges[index], _noises[index], _resolutions[index]);
        }

        #endregion Protected Methods

        public IEnumerable<PossibleProteoform> Deconvolute(double tol, int minEvidence)
        {
            var isotopicPeakGroups = GetIsotopicPeakGroups(this, tol);

            return GetPossibleProteoforms(isotopicPeakGroups, tol);
        }

        private static double ComputeSimilarityScore(IsotopicPeakGroup ok, IsotopicPeakGroup ok2, double tol)
        {
            double score = 0;
            if (ok.charge == ok2.charge)
                return score;
            foreach (var ye in ok.peakList)
                foreach (var ye2 in ok2.peakList)
                    if (ye.X.ToMass(ok.charge) < ye2.X.ToMass(ok2.charge) + tol && ye.X.ToMass(ok.charge) > ye2.X.ToMass(ok2.charge) - tol)
                        score++;
            return score;
        }

        private static IEnumerable<PossibleProteoform> GetPossibleProteoforms(Dictionary<int, List<IsotopicPeakGroup>> isotopicPeakGroups, double tol)
        {
            List<IsotopicPeakGroup> isotopePeakGropusList = new List<IsotopicPeakGroup>();
            foreach (var kvp in isotopicPeakGroups)
                isotopePeakGropusList.AddRange(kvp.Value);

            Dictionary<Tuple<int, int>, double> Similarity = new Dictionary<Tuple<int, int>, double>();
            for (int i = 0; i < isotopePeakGropusList.Count; i++)
            {
                for (int j = i + 1; j < isotopePeakGropusList.Count; j++)
                {
                    Similarity.Add(new Tuple<int, int>(i, j), ComputeSimilarityScore(isotopePeakGropusList[i], isotopePeakGropusList[j], tol));
                }
            }

            PossibleProteoform[,] proteoforms = new PossibleProteoform[isotopePeakGropusList.Count, isotopePeakGropusList.Count];
            HashSet<PossibleProteoform> proteoformsHashSet = new HashSet<PossibleProteoform>();
            foreach (var ye in Similarity.OrderByDescending(b => b.Value).Where(b => b.Value > 0))
            {
                PossibleProteoform a = null;
                for (int i = 0; i < isotopePeakGropusList.Count; i++)
                    if (proteoforms[ye.Key.Item1, i] != null)
                        a = proteoforms[ye.Key.Item1, i];
                PossibleProteoform b = null;
                for (int i = 0; i < isotopePeakGropusList.Count; i++)
                    if (proteoforms[i, ye.Key.Item2] != null)
                        b = proteoforms[i, ye.Key.Item2];

                if (a == null && b == null)
                {
                    PossibleProteoform ok = new PossibleProteoform(isotopePeakGropusList[ye.Key.Item1], isotopePeakGropusList[ye.Key.Item2], ye.Value);
                    proteoforms[ye.Key.Item1, ye.Key.Item2] = ok;
                    proteoforms[ye.Key.Item2, ye.Key.Item1] = ok;
                    proteoformsHashSet.Add(ok);
                    //Console.WriteLine("  new PossibleProteoform ");
                }
                else if (a != null && b == null)
                {
                    proteoforms[ye.Key.Item1, ye.Key.Item2] = a;
                    proteoforms[ye.Key.Item2, ye.Key.Item1] = a;
                    //Console.WriteLine("  adding to existing proteoform: " + a);
                    a.Add(isotopePeakGropusList[ye.Key.Item1], isotopePeakGropusList[ye.Key.Item2], ye.Value);
                    //Console.WriteLine("  added to existing proteoform: " + a);
                }
                else if (a == null && b != null)
                {
                    proteoforms[ye.Key.Item1, ye.Key.Item2] = b;
                    proteoforms[ye.Key.Item2, ye.Key.Item1] = b;
                    //Console.WriteLine("  adding to existing proteoform: " + b);
                    b.Add(isotopePeakGropusList[ye.Key.Item1], isotopePeakGropusList[ye.Key.Item2], ye.Value);
                    //Console.WriteLine("  added to existing proteoform: " + b);
                }
                else if (a != null && b != null)
                {
                    if (a == b)
                    {
                        //Console.WriteLine("  symmetric adding to existing proteoform: " + a);
                        a.Add(isotopePeakGropusList[ye.Key.Item1], isotopePeakGropusList[ye.Key.Item2], ye.Value);
                        //Console.WriteLine("  symmetric adding to existing proteoform: " + a);
                    }
                }
            }



            return proteoformsHashSet;
        }

        private static Dictionary<int, List<IsotopicPeakGroup>> GetIsotopicPeakGroups(ThermoSpectrum mzSpectrum, double tol)
        {
            List<IsotopicPeakGroup> candidatePeakCollections = new List<IsotopicPeakGroup>();
            Dictionary<int, List<IsotopicPeakGroup>> isotopicPeakGroupsByCharge = new Dictionary<int, List<IsotopicPeakGroup>>();
            foreach (ThermoMzPeak peak in mzSpectrum)
            {
                if (peak.Charge > 0)
                {
                    bool peakAccepted = false;
                    if (!isotopicPeakGroupsByCharge.ContainsKey(peak.Charge))
                        isotopicPeakGroupsByCharge.Add(peak.Charge, new List<IsotopicPeakGroup>());
                    foreach (IsotopicPeakGroup ok in isotopicPeakGroupsByCharge[peak.Charge])
                    {
                        if (ok.AttemptToAddNextIsotopicPeak(peak, tol))
                        {
                            peakAccepted = true;
                            break;
                        }
                    }
                    if (!peakAccepted)
                        isotopicPeakGroupsByCharge[peak.Charge].Add(new IsotopicPeakGroup(peak));
                }
            }
            foreach (var ok in isotopicPeakGroupsByCharge)
                isotopicPeakGroupsByCharge[ok.Key].RemoveAll(b => b.Count < 2);
            return isotopicPeakGroupsByCharge;
        }

        public class IsotopicPeakGroup
        {
            public List<ThermoMzPeak> peakList = new List<ThermoMzPeak>();
            public int charge;
            public int Count
            {
                get
                {
                    return peakList.Count;
                }
            }
            public double MostIntenseMass
            {
                get
                {
                    return peakList.OrderByDescending(b => b.Y).First().X.ToMass(charge);
                }

            }
            public override string ToString()
            {
                return "C: " + charge;
            }

            public IsotopicPeakGroup(ThermoMzPeak peak)
            {
                peakList.Add(peak);
                this.charge = peak.Charge;
            }

            internal bool AttemptToAddNextIsotopicPeak(ThermoMzPeak peak, double tol)
            {
                if (peak.Charge == charge &&
                    peak.Mz.ToMass(charge) < (peakList.Last().X.ToMass(charge) + 1) + tol && peak.Mz.ToMass(charge) > (peakList.Last().X.ToMass(charge) + 1) - tol)
                {
                    peakList.Add(peak);
                    return true;
                }
                return false;
            }
        }

        public class PossibleProteoform
        {
            public double EvidenceLevel { get; private set; }

            public HashSet<int> charges = new HashSet<int>();

            public HashSet<IsotopicPeakGroup> isotopicPeakGroups = new HashSet<IsotopicPeakGroup>();

            public double GetMonoisotopicMass()
            {
                double MostIntenseAverage = isotopicPeakGroups.Select(b => b.MostIntenseMass).Average();

                ChemicalFormula chemicalFormula = FindChemicalFormulaWithAverageMass(MostIntenseAverage);

                var ye = IsotopicDistribution.GetDistribution(chemicalFormula, 0.01, 0.001);
                double[] massesArray = ye.Masses.ToArray();
                double[] intensitiesArray = ye.Intensities.ToArray();
                Array.Sort(intensitiesArray, massesArray);
                double shiftFromMostIntenseToMono = massesArray.Last() - chemicalFormula.MonoisotopicMass;

                return MostIntenseMode() - shiftFromMostIntenseToMono;
            }

            private double MostIntenseMode()
            {
                Dictionary<double, List<double>> mostIntenses = new Dictionary<double, List<double>>();
                foreach (var ye in isotopicPeakGroups)
                {
                    bool added = false;
                    foreach (var kvp in mostIntenses)
                    {
                        if (kvp.Key > ye.MostIntenseMass - 0.5 && kvp.Key < ye.MostIntenseMass + 0.5)
                        {
                            kvp.Value.Add(ye.MostIntenseMass);
                            added = true;
                            break;
                        }
                    }
                    if (!added)
                        mostIntenses.Add(ye.MostIntenseMass, new List<double>() { ye.MostIntenseMass });
                }
                int bestCount = 0;
                KeyValuePair<double, List<double>> bestKvp = new KeyValuePair<double, List<double>>();
                foreach (var kvp in mostIntenses)
                {
                    if (kvp.Value.Count > bestCount)
                    {
                        bestCount = kvp.Value.Count;
                        bestKvp = kvp;
                    }
                }
                return bestKvp.Value.Average();
            }


            private ChemicalFormula FindChemicalFormulaWithAverageMass(double mostIntenseAverage)
            {
                const double averageC = 4.9384;
                const double averageH = 7.7583;
                const double averageO = 1.4773;
                const double averageN = 1.3577;
                const double averageS = 0.0417;

                ChemicalFormula chemicalFormula = new ChemicalFormula();

                double factor = (mostIntenseAverage / 98.123);

                chemicalFormula.Add("C", Convert.ToInt32(Math.Round(averageC * factor)));
                chemicalFormula.Add("H", Convert.ToInt32(Math.Round(averageH * factor)));
                chemicalFormula.Add("O", Convert.ToInt32(Math.Round(averageO * factor)));
                chemicalFormula.Add("N", Convert.ToInt32(Math.Round(averageN * factor)));
                chemicalFormula.Add("S", Convert.ToInt32(Math.Round(averageS * factor)));

                do
                {
                    if (chemicalFormula.AverageMass > mostIntenseAverage + 33 && chemicalFormula.ContainsIsotopesOf(PeriodicTable.GetElement("S")))
                        chemicalFormula.Remove(ChemicalFormula.ParseFormula("S"));
                    if (chemicalFormula.AverageMass > mostIntenseAverage + 17 && chemicalFormula.ContainsIsotopesOf(PeriodicTable.GetElement("O")))
                        chemicalFormula.Remove(ChemicalFormula.ParseFormula("N"));
                    if (chemicalFormula.AverageMass > mostIntenseAverage + 15 && chemicalFormula.ContainsIsotopesOf(PeriodicTable.GetElement("N")))
                        chemicalFormula.Remove(ChemicalFormula.ParseFormula("O"));
                    if (chemicalFormula.AverageMass > mostIntenseAverage + 13 && chemicalFormula.ContainsIsotopesOf(PeriodicTable.GetElement("C")))
                        chemicalFormula.Remove(ChemicalFormula.ParseFormula("C"));
                    if (chemicalFormula.AverageMass > mostIntenseAverage + 1.001 && chemicalFormula.ContainsIsotopesOf(PeriodicTable.GetElement("H")))
                        chemicalFormula.Remove(ChemicalFormula.ParseFormula("H"));

                    if (chemicalFormula.AverageMass < mostIntenseAverage - 33)
                        chemicalFormula.Add(ChemicalFormula.ParseFormula("S"));
                    if (chemicalFormula.AverageMass < mostIntenseAverage - 17)
                        chemicalFormula.Add(ChemicalFormula.ParseFormula("N"));
                    if (chemicalFormula.AverageMass < mostIntenseAverage - 15)
                        chemicalFormula.Add(ChemicalFormula.ParseFormula("O"));
                    if (chemicalFormula.AverageMass < mostIntenseAverage - 13)
                        chemicalFormula.Add(ChemicalFormula.ParseFormula("C"));
                    if (chemicalFormula.AverageMass < mostIntenseAverage - 1.001)
                        chemicalFormula.Add(ChemicalFormula.ParseFormula("H"));
                } while ((chemicalFormula.AverageMass > mostIntenseAverage + 1.001) || (chemicalFormula.AverageMass < mostIntenseAverage - 1.001));

                return chemicalFormula;
            }

            public PossibleProteoform(IsotopicPeakGroup isotopicPeakGroup1, IsotopicPeakGroup isotopicPeakGroup2, double score)
            {
                isotopicPeakGroups.Add(isotopicPeakGroup1);
                isotopicPeakGroups.Add(isotopicPeakGroup2);
                charges.Add(isotopicPeakGroup1.charge);
                charges.Add(isotopicPeakGroup2.charge);
                EvidenceLevel = score;
            }

            public override string ToString()
            {
                return "MonoisotopicMass: " + GetMonoisotopicMass() + " charges: " + string.Join(",", charges.OrderBy(b => b)) + " evidence level: " + EvidenceLevel;
            }

            internal void Add(IsotopicPeakGroup isotopicPeakGroup1, IsotopicPeakGroup isotopicPeakGroup2, double score)
            {
                isotopicPeakGroups.Add(isotopicPeakGroup1);
                isotopicPeakGroups.Add(isotopicPeakGroup2);
                charges.Add(isotopicPeakGroup1.charge);
                charges.Add(isotopicPeakGroup2.charge);
                EvidenceLevel += score;
            }
        }



    }
}