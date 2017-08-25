﻿// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (MzSpectrum.cs) is part of MassSpectrometry.
//
// MassSpectrometry is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MassSpectrometry is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with MassSpectrometry. If not, see <http://www.gnu.org/licenses/>.

using Chemistry;
using MzLibUtil;
using Spectra;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace MassSpectrometry
{
    public abstract class MzSpectrum<TPeak> : Spectrum<TPeak>, IMzSpectrum<TPeak>
        where TPeak : IMzPeak
    {
        #region Private Fields

        private static readonly double[] mms = new double[] { 1.0029, 2.0052, 3.0077, 4.01, 5.012, 6.0139, 7.0154, 8.0164 };

        private static readonly List<Tuple<double, List<double>>> intensityFractions = new List<Tuple<double, List<double>>>();

        #endregion Private Fields

        #region Public Constructors

        static MzSpectrum()
        {
            intensityFractions.Add(new Tuple<double, List<double>>(155, new List<double> { 0.915094568, 0.07782302, 0.006528797, 0.000289506 }));
            intensityFractions.Add(new Tuple<double, List<double>>(226, new List<double> { 0.88015657, 0.107467263, 0.011417303, 0.000730494 }));
            intensityFractions.Add(new Tuple<double, List<double>>(310, new List<double> { 0.837398069, 0.142430845, 0.01821746, 0.001683771 }));
            intensityFractions.Add(new Tuple<double, List<double>>(437, new List<double> { 0.777595132, 0.186958768, 0.031114269, 0.003704342, 0.000220493 }));
            intensityFractions.Add(new Tuple<double, List<double>>(620, new List<double> { 0.701235526, 0.238542629, 0.050903269, 0.008082801, 0.000985192 }));
            intensityFractions.Add(new Tuple<double, List<double>>(888, new List<double> { 0.602453248, 0.291899044, 0.084076553, 0.01790019, 0.002916629, 0.000410371 }));
            intensityFractions.Add(new Tuple<double, List<double>>(1243, new List<double> { 0.492328432, 0.333344333, 0.128351944, 0.035959923, 0.008063481, 0.001433271, 0.000195251 }));
            intensityFractions.Add(new Tuple<double, List<double>>(1797, new List<double> { 0.348495022, 0.336686099, 0.193731423, 0.082270917, 0.028068866, 0.008052644, 0.001907311, 0.000372359, 4.52281E-05 }));
            intensityFractions.Add(new Tuple<double, List<double>>(2515, new List<double> { 0.229964408, 0.313975523, 0.238643189, 0.130654102, 0.056881604, 0.020732138, 0.006490044, 0.001706308, 0.000373761, 4.55951E-05 }));
            intensityFractions.Add(new Tuple<double, List<double>>(3532, new List<double> { 0.12863395, 0.247015676, 0.254100853, 0.184302695, 0.104989402, 0.049731171, 0.020279668, 0.007267861, 0.002300006, 0.000619357, 9.64322E-05 }));
            intensityFractions.Add(new Tuple<double, List<double>>(5019, new List<double> { 0.053526677, 0.145402081, 0.208920636, 0.209809764, 0.164605485, 0.107024765, 0.059770563, 0.029447041, 0.012957473, 0.005127018, 0.001845335, 0.000572486, 0.000115904 }));
        }

        #endregion Public Constructors

        #region Protected Constructors

        protected MzSpectrum(double[,] mzintensities) : base(mzintensities)
        {
        }

        protected MzSpectrum(double[] mz, double[] intensities, bool shouldCopy) : base(mz, intensities, shouldCopy)
        {
        }

        #endregion Protected Constructors

        #region Public Properties

        new public MzRange Range
        {
            get
            {
                return new MzRange(FirstX, LastX);
            }
        }

        #endregion Public Properties

        #region Public Methods

        public static byte[] Get64Bitarray(IEnumerable<double> yArray)
        {
            var mem = new MemoryStream();
            foreach (var okk in yArray)
            {
                byte[] ok = BitConverter.GetBytes(okk);
                mem.Write(ok, 0, ok.Length);
            }
            mem.Position = 0;
            return mem.ToArray();
        }

        public byte[] Get64BitYarray()
        {
            return Get64Bitarray(YArray);
        }

        public byte[] Get64BitXarray()
        {
            return Get64Bitarray(XArray);
        }

        public void ReplaceXbyApplyingFunction(Func<IMzPeak, double> convertor)
        {
            for (int i = 0; i < Size; i++)
                XArray[i] = convertor(this[i]);
            peakWithHighestY = default(TPeak);
            peakList = new TPeak[Size];
        }

        public override string ToString()
        {
            return string.Format("{0} (Peaks {1})", Range, Size);
        }

        public IEnumerable<Tuple<List<IMzPeak>, int>> Deconvolute(MzRange theRange, int maxAssumedChargeState, Tolerance massTolerance, double intensityRatio)
        {
            var isolatedMassesAndCharges = new List<Tuple<List<IMzPeak>, int>>();

            foreach (var peak in Extract(theRange))
            {
                // Always assume the current peak is a monoisotopic peak!

                List<IMzPeak> bestListOfPeaks = new List<IMzPeak>();
                int bestChargeState = 1;
                for (int chargeState = 1; chargeState <= maxAssumedChargeState; chargeState++)
                {
                    var listOfPeaksForThisChargeState = new List<IMzPeak> { peak };
                    var mMass = peak.Mz.ToMass(chargeState);
                    for (int mm = 1; mm <= mms.Length; mm++)
                    {
                        double diffToNextMmPeak = mms[mm - 1];
                        double theorMass = mMass + diffToNextMmPeak;
                        var closestpeak = GetClosestPeak(theorMass.ToMz(chargeState));
                        if (massTolerance.Within(closestpeak.Mz.ToMass(chargeState), theorMass) && SatisfiesRatios(mMass, mm, peak, closestpeak, intensityRatio))
                        {
                            // Found a match to an isotope peak for this charge state!
                            listOfPeaksForThisChargeState.Add(closestpeak);
                        }
                        else
                            break;
                    }
                    if (listOfPeaksForThisChargeState.Count >= bestListOfPeaks.Count)
                    {
                        bestListOfPeaks = listOfPeaksForThisChargeState;
                        bestChargeState = chargeState;
                    }
                }
                if (bestListOfPeaks.Count >= 2)
                    isolatedMassesAndCharges.Add(new Tuple<List<IMzPeak>, int>(bestListOfPeaks, bestChargeState));
            }

            List<double> seen = new List<double>();
            while (isolatedMassesAndCharges.Any())
            {
                // Pick longest
                var longest = isolatedMassesAndCharges.OrderByDescending(b => b.Item1.Count).First();
                yield return longest;
                isolatedMassesAndCharges.Remove(longest);
                isolatedMassesAndCharges.RemoveAll(b => b.Item1.Intersect(longest.Item1).Any());
            }
        }

        #endregion Public Methods

        #region Private Methods

        private bool SatisfiesRatios(double mMass, int mm, IMzPeak ye, IMzPeak closestpeak, double intensityRatio)
        {
            double bestDiff = double.MaxValue;
            List<double> bestFracList = null;
            for (int i = 0; i < intensityFractions.Count; i++)
            {
                var diff = Math.Abs(mMass - intensityFractions[i].Item1);
                if (diff < bestDiff)
                {
                    bestDiff = diff;
                    bestFracList = intensityFractions[i].Item2;
                }
            }
            if (bestFracList == null || bestFracList.Count <= mm)
                return false;

            var theMM = bestFracList[0];
            var theCompared = bestFracList[mm];

            var comparedShouldBe = ye.Intensity / theMM * theCompared;

            if (closestpeak.Intensity < comparedShouldBe / intensityRatio || closestpeak.Intensity > comparedShouldBe * intensityRatio)
                return false;

            return true;
        }

        #endregion Private Methods
    }
}