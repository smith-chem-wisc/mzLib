using Chemistry;
using Easy.Common.Extensions;
using MassSpectrometry;
using MathNet.Numerics;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Documents;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class TestShortreedConv
    {

        [Test]
        public static void Bubba()
        {
            Dictionary<int, string> hashToPeptide = new Dictionary<int, string>();
            Dictionary<int, List<int>> peakToListHash = new Dictionary<int, List<int>>();

            List<Protein> proteins = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, "DeconTests", @"hela_snip_for_unitTest.fasta"), true, DecoyType.None, false, out var a,
                ProteinDbLoader.EnsemblAccessionRegex, ProteinDbLoader.EnsemblFullNameRegex, ProteinDbLoader.EnsemblAccessionRegex, ProteinDbLoader.EnsemblGeneNameRegex, null);

            foreach (Protein protein in proteins)
            {
                foreach (PeptideWithSetModifications pwsm in protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()))
                {
                    int peptideHash = pwsm.GetHashCode();
                    if (!hashToPeptide.ContainsKey(peptideHash))
                    {
                        hashToPeptide.Add(peptideHash, pwsm.ToString());
                    }

                    foreach (int mz in MzValuesSignificant(pwsm))
                    {
                        if (peakToListHash.ContainsKey(mz))
                        {
                            peakToListHash[mz].Add(peptideHash);
                        }
                        else
                        {
                            peakToListHash.Add(mz, new List<int> { peptideHash });
                        }
                    }

                }
            }

            foreach (KeyValuePair<int, List<int>> item in peakToListHash)
            {
                peakToListHash[item.Key] = item.Value.Distinct().ToList();
            }

            var reader = MsDataFileReader.GetDataFile(Path.Combine(TestContext.CurrentContext.TestDirectory,
                    "DeconTests", "TaGe_SA_HeLa_04_subset_longestSeq.mzML"));
            reader.LoadAllStaticData();

            var ms1scans = reader.Scans.Where(s=>s.MsnOrder == 1).ToList();

            List<string> myOUt = new List<string>();

            foreach (MsDataScan spectrum in ms1scans)
            {
                double[] intensities = spectrum.MassSpectrum.YArray;
                double[] mzs = spectrum.MassSpectrum.XArray;

                Array.Sort(intensities, mzs);

                Array.Reverse(mzs,0,mzs.Length);

                int max = Math.Min(400,mzs.Length);

                List<int> mzsToLookUp = new List<int>();

                for (int i = 0; i < max; i++)
                {
                    mzsToLookUp.Add((int)(mzs[i] * 100).Round(0));
                }

                List<int> hashes = new List<int>();

                foreach (var mz in mzsToLookUp)
                {
                    if (peakToListHash.ContainsKey(mz))
                    {
                        hashes.AddRange(peakToListHash[mz]);
                    }
                }

                var most = hashes.GroupBy(i => i).OrderByDescending(grp => grp.Count()).Select(grp => grp.Key).ToList();

                for (int i = 0; i < 10; i++)
                {
                    if (hashToPeptide.ContainsKey(most[i]))
                    {
                        myOUt.Add(spectrum.OneBasedScanNumber + "\t" + hashToPeptide[most[i]]);
                    }
                }
  
            }
            File.WriteAllLines(@"C:\Users\Michael Shortreed\Downloads\try.txt", myOUt.ToArray());
            Assert.IsTrue(false);
        }

        public static List<int> MzValuesSignificant(PeptideWithSetModifications pwsm)
        {
            var result = new List<int>();

            ChemicalFormula formula = pwsm.FullChemicalFormula;

            var isotopicDistribution = IsotopicDistribution.GetDistribution(formula, 0.125, 0.0001);

            double[] masses = isotopicDistribution.Masses.ToArray();
            double[] intensities = isotopicDistribution.Intensities.ToArray();

            Array.Sort(intensities, masses);
            Array.Reverse(masses,0,masses.Length);

            int max = Math.Min(50,masses.Length);

            List<double> acceptableMasses = masses[0..max].ToList();

            for (int i = 2; i < 5; i++)
            {
                foreach (var mass in acceptableMasses)
                {
                    result.Add((int)(mass.ToMz(i) * 100).Round(0));
                }
            }

            return result;
        }
    }
}
