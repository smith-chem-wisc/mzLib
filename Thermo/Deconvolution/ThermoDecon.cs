using Chemistry;
using System;
using System.Collections.Generic;
using System.Linq;

namespace IO.Thermo.Deconvolution
{
    internal class ThermoDecon
    {

        #region Internal Methods

        internal static IEnumerable<PossibleProteoform> GetPossibleProteoforms(Dictionary<int, List<IsotopicPeakGroup>> isotopicPeakGroups, double tol)
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

        internal static Dictionary<int, List<IsotopicPeakGroup>> GetIsotopicPeakGroups(ThermoSpectrum mzSpectrum, double tol)
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

        #endregion Internal Methods

        #region Private Methods

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

        #endregion Private Methods

    }
}