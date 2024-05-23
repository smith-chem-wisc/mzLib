// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work Copyright 2016 Stefan Solntsev
//
// This file (TestSpectra.cs) is part of MassSpectrometry.Tests.
//
// MassSpectrometry.Tests is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MassSpectrometry.Tests is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with MassSpectrometry.Tests. If not, see <http://www.gnu.org/licenses/>.

using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics.PSM;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestSpectra
    {
        private MzSpectrum _mzSpectrumA;
        private static Stopwatch Stopwatch { get; set; }

        [SetUp]
        public static void Setuppp()
        {
            Stopwatch = new Stopwatch();
            Stopwatch.Start();
        }

        [TearDown]
        public static void TearDown()
        {
            Console.WriteLine($"Analysis time: {Stopwatch.Elapsed.Hours}h {Stopwatch.Elapsed.Minutes}m {Stopwatch.Elapsed.Seconds}s");
        }

        [SetUp]
        public void Setup()
        {
            double[] mz = { 328.73795, 329.23935, 447.73849, 448.23987, 482.23792, 482.57089, 482.90393, 500.95358, 501.28732, 501.62131, 611.99377, 612.32806, 612.66187, 722.85217, 723.35345 };
            double[] intensities = { 81007096.0, 28604418.0, 78353512.0, 39291696.0, 122781408.0, 94147520.0, 44238040.0, 71198680.0, 54184096.0, 21975364.0, 44514172.0, 43061628.0, 23599424.0, 56022696.0, 41019144.0 };

            _mzSpectrumA = new MzSpectrum(mz, intensities, false);
        }



        [Test]
        public static void MoreJunk()
        {
            List<string> quantResults = File.ReadAllLines(@"C:\Users\Michael Shortreed\Downloads\4-26-24 AD Brain DHAA Analysis\AllQuantifiedPeptidesMarkovich.tsv").ToList();
            List<string> fullSequences = new List<string>();
            foreach (var line in quantResults)
            {
                string[] fields = line.Split('\t');
                fullSequences.Add(fields[0]);
            }

            string psmFilePath =
                @"C:\Users\Michael Shortreed\Downloads\4-26-24 AD Brain DHAA Analysis\AllPeptidesMarkovich.psmtsv";
            List<PsmFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFilePath, out var warnings);

            List<PsmFromTsv> quantPeptidesFoundInAllPeptides = new List<PsmFromTsv>();
            foreach (var psm in parsedPsms)
            {
                bool containsAny = fullSequences.Any(s => psm.FullSequence.Contains(s));
                if (containsAny)
                {
                    fullSequences.Remove(psm.FullSequence);
                    quantPeptidesFoundInAllPeptides.Add(psm);
                }

                if (fullSequences.Count == 0)
                {
                    break;
                }
            }

            List<string> myOut = new List<string>();
            //myOut.Add(parsedPsms[0].ToString());
            foreach (var psm in quantPeptidesFoundInAllPeptides)
            {

                myOut.Add(psm.ToString());

            }

            File.WriteAllLines(@"C:\Users\Michael Shortreed\Downloads\4-26-24 AD Brain DHAA Analysis\psmsFoundInQuantMarkovich.txt", myOut);

        }
        [Test]
        public static void Junk3()
        {

            string psmFilePath =
                @"C:\Users\Michael Shortreed\Downloads\4-26-24 AD Brain DHAA Analysis\AllPeptidesMarkovich.psmtsv";
            List<PsmFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFilePath, out var warnings);
            parsedPsms = parsedPsms.Where(p => p.QValue < .01).ToList();
            parsedPsms = parsedPsms.Where(p => p.DecoyContamTarget.Contains("T")).ToList();

            List<string> interestingMods = new List<string> { "Y[Common Biological:Phosphorylation on Y]", "S[Common Biological:Phosphorylation on S]", "T[Common Biological:Phosphorylation on T]" };

            Dictionary<(string, int), List<PsmFromTsv>> genePositionPsm = new Dictionary<(string, int), List<PsmFromTsv>>();

            foreach (var psm in parsedPsms)
            {
                List<string> foundMods = interestingMods.Where(s => psm.FullSequence.Contains(s)).ToList();
                if (foundMods.Any())
                {
                    string sequence = psm.FullSequence;
                    foreach (var mod in foundMods)
                    {
                        string firstCharacter = mod.Substring(0, 1);
                        if (firstCharacter == "Y" || firstCharacter == "S" || firstCharacter == "T")
                        {
                            sequence = sequence.Replace(mod, firstCharacter.ToLowerInvariant());
                        }
                    }
                    sequence = Regex.Replace(sequence, "\\[(.*?)\\]", "");
                    List<int> lowercaseLetterPositions = new List<int>();

                    for (int i = 0; i < sequence.Length; i++)
                    {
                        if (char.IsLower(sequence[i]))
                        {
                            lowercaseLetterPositions.Add(i);
                        }
                    }
                    string firstAndLastAminoAcidPositionInProtein = psm.StartAndEndResiduesInProtein.Split('|')[0];
                    firstAndLastAminoAcidPositionInProtein = firstAndLastAminoAcidPositionInProtein.Substring(1, firstAndLastAminoAcidPositionInProtein.Length - 2);
                    firstAndLastAminoAcidPositionInProtein = firstAndLastAminoAcidPositionInProtein.Replace(" to ", "\t");
                    int[] startEnd = firstAndLastAminoAcidPositionInProtein.Split('\t').Select(int.Parse).ToArray();
                    lowercaseLetterPositions = lowercaseLetterPositions.Select(s => s + startEnd[0]).ToList();

                    string allGenesInPsm = psm.GeneName;
                    string[] genes = allGenesInPsm.Split('|');
                    string firstGene = genes[0].Split(':')[1];

                    foreach (int position in lowercaseLetterPositions)
                    {
                        if (genePositionPsm.ContainsKey((firstGene, position)))
                        {
                            genePositionPsm[(firstGene, position)].Add(psm);

                        }
                        else
                        {
                            genePositionPsm.Add((firstGene, position), new List<PsmFromTsv> { psm });
                        }
                    }
                }
            }

            List<string> myOut = new List<string>();
            myOut.Add("position" + "\t" + "protein accession" + "\t" + "gene" + "\t" + "full sequence");

            foreach (var kvp in genePositionPsm)
            {
                if (kvp.Value.Count > 1)
                {
                    foreach (var psm in kvp.Value)
                    {
                        myOut.Add(kvp.Key.Item2 + "\t" + psm.ProteinAccession + "\t" + kvp.Key.Item1 + "\t" + psm.FullSequence);
                    }
                }
            }

            File.WriteAllLines(@"C:\Users\Michael Shortreed\Downloads\4-26-24 AD Brain DHAA Analysis\phosphoPeptidesMarkovich.txt", myOut);

        }

        [Test]
        public static void Junk4()
        {
            Dictionary<(string,string),string> baseSequencefullSequencQvaluePepQvalueforPSMs = new Dictionary<(string, string), string>();

            using (StreamReader sr = new StreamReader(@"C:\Users\Michael Shortreed\Downloads\4-26-24 AD Brain DHAA Analysis\AllPSMsMarkovich.psmtsv"))
            {
                bool continueReading = true;
                string line;
                while ((line = sr.ReadLine()) != null && continueReading)
                {
                    string[] fields = line.Split('\t');
                    string baseSequence = fields[12];
                    string fullSequence = fields[13];
                    string targetDecoyContam = fields[38];
                    string qValue = fields[50];
                    string pepQValue = fields[55];
                    if (targetDecoyContam == "T" && double.TryParse(qValue, out double qValueDouble))
                    {
                        if (qValueDouble < 0.01 &&!baseSequencefullSequencQvaluePepQvalueforPSMs.ContainsKey((baseSequence,fullSequence)))
                        {
                            baseSequencefullSequencQvaluePepQvalueforPSMs.Add((baseSequence, fullSequence), qValue + "\t" + pepQValue);
                        }

                    }

                    if (double.TryParse(qValue, out double qValueDouble2))
                    {
                        if (qValueDouble2 > 0.01)
                        {
                            continueReading = false;
                        }
                    }
                }
            }

            string psmFilePath =
                @"C:\Users\Michael Shortreed\Downloads\4-26-24 AD Brain DHAA Analysis\AllPeptidesMarkovich.psmtsv";
            List<PsmFromTsv> parsedPeptides = SpectrumMatchTsvReader.ReadPsmTsv(psmFilePath, out var warnings);
            parsedPeptides = parsedPeptides.Where(p => p.DecoyContamTarget.Contains("T")).ToList();

            List<string> interestingMods = new List<string> { "S[Common Biological:Phosphorylation on S]", "T[Common Biological:Phosphorylation on T]", "S[Less Common:Dehydroalanine on S]", "C[Less Common:Dehydroalanine on C]",
                "S[Custom:Homocys on S]", "C[Custom:Homocys on C]", "T[Custom:Homocys on T]", "T[Less Common:Dehydrobutyrine on T]", "S[Custom:DTT on S]", "C[Custom:DTT on C]", "C[Custom:DTT on T]", "T[Custom:Glutathione on T]", 
                "S[Custom:Glutathione on S]", "C[Custom:Glutathione on C]", "S[Custom:TCEP on S]", "T[Custom:TCEP on T]", "C[Custom:TCEP on C]"  };

            Dictionary<(string, int), List<PsmFromTsv>> genePositionPeptides = new Dictionary<(string, int), List<PsmFromTsv>>();
            Dictionary<(string, int), List<string>> genePositionMod = new Dictionary<(string, int), List<string>>();
            Dictionary<string,List<PsmFromTsv>> noInterestingModsPeptideValues = new Dictionary<string, List<PsmFromTsv>>();

            foreach (var peptide in parsedPeptides)
            {
                if (baseSequencefullSequencQvaluePepQvalueforPSMs.ContainsKey((peptide.BaseSeq,peptide.FullSequence)))
                {
                    List<string> foundMods = interestingMods.Where(s => peptide.FullSequence.Contains(s)).ToList();
                    if (foundMods.Any())
                    {
                        foreach (var mod in foundMods)
                        {
                            string sequence = peptide.FullSequence;
                            string firstCharacter = mod.Substring(0, 1);
                            if (firstCharacter == "Y" || firstCharacter == "S" || firstCharacter == "T" ||
                                firstCharacter == "C")
                            {
                                sequence = sequence.Replace(mod, firstCharacter.ToLowerInvariant());
                            }

                            //eliminate the remaning mods
                            while (sequence.Contains("[") && sequence.Contains("]"))
                            {
                                int firstOpenBracket = sequence.IndexOf('[');
                                int firstCloseBracket = sequence.IndexOf(']', firstOpenBracket);
                                if (firstCloseBracket != -1)
                                {
                                    sequence = sequence.Remove(firstOpenBracket, firstCloseBracket - firstOpenBracket + 1);
                                }
                                else
                                {
                                    break;
                                }
                            }

                            List<int> lowercaseLetterPositions = new List<int>();

                            for (int i = 0; i < sequence.Length; i++)
                            {
                                if (char.IsLower(sequence[i]))
                                {
                                    lowercaseLetterPositions.Add(i);
                                }
                            }

                            string firstAndLastAminoAcidPositionInProtein =
                                peptide.StartAndEndResiduesInProtein.Split('|')[0];
                            firstAndLastAminoAcidPositionInProtein =
                                firstAndLastAminoAcidPositionInProtein.Substring(1,
                                    firstAndLastAminoAcidPositionInProtein.Length - 2);
                            firstAndLastAminoAcidPositionInProtein =
                                firstAndLastAminoAcidPositionInProtein.Replace(" to ", "\t");
                            int[] startEnd = firstAndLastAminoAcidPositionInProtein.Split('\t').Select(int.Parse)
                                .ToArray();
                            lowercaseLetterPositions = lowercaseLetterPositions.Select(s => s + startEnd[0]).ToList();

                            string allGenesInPsm = peptide.GeneName;
                            string firstGene = "";
                            if (allGenesInPsm.Contains("|"))
                            {
                                string[] genes = allGenesInPsm.Split('|');

                                if (genes[0].Contains(":"))
                                {
                                    firstGene = genes[0].Split(':')[1];

                                    foreach (int position in lowercaseLetterPositions)
                                    {
                                        if (genePositionPeptides.ContainsKey((firstGene, position)))
                                        {
                                            genePositionPeptides[(firstGene, position)].Add(peptide);
                                            genePositionMod[(firstGene, position)]
                                                .Add(mod + "\t" + peptide.QValue + "\t" + peptide.PEP_QValue);

                                        }
                                        else
                                        {
                                            genePositionPeptides.Add((firstGene, position), new List<PsmFromTsv> { peptide });
                                            genePositionMod.Add((firstGene, position),
                                                new List<string> { mod + "\t" + peptide.QValue + "\t" + peptide.PEP_QValue });
                                        }
                                    }
                                }
                                else
                                {
                                    firstGene = genes[0];
                                    if (genes[0].Contains(":"))
                                    {
                                        firstGene = genes[0].Split(':')[1];
                                    }
                                    foreach (int position in lowercaseLetterPositions)
                                    {
                                        if (genePositionPeptides.ContainsKey((firstGene, position)))
                                        {
                                            genePositionPeptides[(firstGene, position)].Add(peptide);
                                            genePositionMod[(firstGene, position)]
                                                .Add(mod + "\t" + peptide.QValue + "\t" + peptide.PEP_QValue);

                                        }
                                        else
                                        {
                                            genePositionPeptides.Add((firstGene, position), new List<PsmFromTsv> { peptide });
                                            genePositionMod.Add((firstGene, position),
                                                new List<string> { mod + "\t" + peptide.QValue + "\t" + peptide.PEP_QValue });
                                        }
                                    }
                                }
                            }
                            else
                            {
                                firstGene = peptide.GeneName;
                                if (peptide.GeneName.Contains(":"))
                                {
                                    firstGene = peptide.GeneName.Split(':')[1];
                                }

                                foreach (int position in lowercaseLetterPositions)
                                {
                                    if (genePositionPeptides.ContainsKey((firstGene, position)))
                                    {
                                        genePositionPeptides[(firstGene, position)].Add(peptide);
                                        genePositionMod[(firstGene, position)]
                                            .Add(mod + "\t" + peptide.QValue + "\t" + peptide.PEP_QValue);

                                    }
                                    else
                                    {
                                        genePositionPeptides.Add((firstGene, position), new List<PsmFromTsv> { peptide });
                                        genePositionMod.Add((firstGene, position),
                                            new List<string> { mod + "\t" + peptide.QValue + "\t" + peptide.PEP_QValue });
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        if (noInterestingModsPeptideValues.ContainsKey(peptide.BaseSeq))
                        {
                            noInterestingModsPeptideValues[peptide.BaseSeq].Add(peptide);
                        }
                        else
                        {
                            noInterestingModsPeptideValues.Add(peptide.BaseSeq, new List<PsmFromTsv> { peptide });
                        }
                        
                    }
                }
            }

            List<string> myOut = new List<string>();
            myOut.Add("position" + "\t" + "protein accession" + "\t" + "gene" + "\t" + "base sequence" + "\t" + "full sequence" + "\t" + "modification" + "\t" + "Peptide Q-value" + "\t" + "Peptide PEP Q-Value" + "\t" + "PSM Q-value" + "\t" + "PSM PEP Q-Value");

            foreach (var kvp in genePositionPeptides)
            {
                if (kvp.Value.Count > 1)
                {
                    foreach (var psm in kvp.Value)
                    {
                        int index = kvp.Value.IndexOf(psm);
                        if(baseSequencefullSequencQvaluePepQvalueforPSMs.ContainsKey((psm.BaseSeq, psm.FullSequence)))
                        {
                            myOut.Add(kvp.Key.Item2 + "\t" + psm.ProteinAccession + "\t" + kvp.Key.Item1 + "\t" + psm.BaseSeq + "\t" + psm.FullSequence + "\t" + genePositionMod[kvp.Key][index] + "\t" + baseSequencefullSequencQvaluePepQvalueforPSMs[(psm.BaseSeq, psm.FullSequence)]);
                            if (noInterestingModsPeptideValues.ContainsKey(psm.BaseSeq))
                            {
                                foreach (var psm2 in noInterestingModsPeptideValues[psm.BaseSeq])
                                {
                                    myOut.Add(kvp.Key.Item2 + "\t" + psm2.ProteinAccession + "\t" + kvp.Key.Item1 + "\t" + psm2.BaseSeq + "\t" + psm2.FullSequence + "\t" + "No Modifications of Interest" + "\t" + psm2.QValue + "\t" + psm2.PEP_QValue + "\t" + baseSequencefullSequencQvaluePepQvalueforPSMs[(psm2.BaseSeq, psm2.FullSequence)]);
                                }
                                noInterestingModsPeptideValues.Remove(psm.BaseSeq);
                            }
                        }
                        else
                        {
                            myOut.Add(kvp.Key.Item2 + "\t" + psm.ProteinAccession + "\t" + kvp.Key.Item1 + "\t" + psm.BaseSeq + "\t" + psm.FullSequence + "\t" + genePositionMod[kvp.Key][index]);
                            if (noInterestingModsPeptideValues.ContainsKey(psm.BaseSeq))
                            {
                                foreach (var psm2 in noInterestingModsPeptideValues[psm.BaseSeq])
                                {
                                    string psmQvalues = baseSequencefullSequencQvaluePepQvalueforPSMs[(psm2.BaseSeq, psm2.FullSequence)];
                                    myOut.Add(kvp.Key.Item2 + "\t" + psm2.ProteinAccession + "\t" + kvp.Key.Item1 + "\t" + psm2.BaseSeq + "\t" + psm2.FullSequence + "\t" + "No Modifications of Interest" + "\t" + psm2.QValue + "\t" + psm2.PEP_QValue + "\t" + psmQvalues);
                                
                                }
                                noInterestingModsPeptideValues.Remove(psm.BaseSeq);
                            }
                        }
                        
                    }
                }
            }

            File.WriteAllLines(@"C:\Users\Michael Shortreed\Downloads\4-26-24 AD Brain DHAA Analysis\modPeptidesMarkovichWithPosition.txt", myOut);

        }



        [Test]
        public static void Junk5()
        {
            List<string> interestingMods = new List<string> { "S[Common Biological:Phosphorylation on S]", "T[Common Biological:Phosphorylation on T]", "S[Less Common:Dehydroalanine on S]", "C[Less Common:Dehydroalanine on C]",
                "S[Custom:Homocys on S]", "C[Custom:Homocys on C]", "T[Custom:Homocys on T]", "T[Less Common:Dehydrobutyrine on T]", "S[Custom:DTT on S]", "C[Custom:DTT on C]", "C[Custom:DTT on T]", "T[Custom:Glutathione on T]",
                "S[Custom:Glutathione on S]", "C[Custom:Glutathione on C]", "S[Custom:TCEP on S]", "T[Custom:TCEP on T]", "C[Custom:TCEP on C]"  };

            //string is position accession gene
            //dictionary key is mod and int is count
            Dictionary<string, Dictionary<string, int>> bubba = new Dictionary<string, Dictionary<string, int>>();
            using (StreamReader sr =
                   new StreamReader(
                       @"C:\Users\Michael Shortreed\Downloads\4-26-24 AD Brain DHAA Analysis\modPeptidesMarkovichWithPosition.txt"))
            {
                bool continueReading = true;
                string line;
                
                while ((line = sr.ReadLine()) != null && continueReading)
                {
                    
                    string[] fields = line.Split('\t');
                    string myKey = fields[0] + "\t" + fields[1] + "\t" + fields[2];

                    bool goodLine = (double.TryParse(fields[5], out double qvalue) && qvalue < 0.01 );

                    if (goodLine)
                    {
                        if (bubba.ContainsKey(myKey))
                        {
                            if (bubba[myKey].ContainsKey(fields[4]))
                            {
                                bubba[myKey][fields[4]]++;
                            }
                            else
                            {
                                bubba[myKey].Add(fields[4], 1);
                            }
                        }
                        else
                        {
                            bubba.Add(myKey, new Dictionary<string, int> { { fields[4], 1 } });
                        }
                    }
                }
            }

            List<string> myOut = new List<string>();
            myOut.Add("position" + "\t" + "protein accession" + "\t" + "gene" + "\t" + String.Join('\t',interestingMods));

            foreach (var kvp in bubba)
            {
                string myLine = kvp.Key +"\t";
                foreach (var mod in interestingMods)
                {
                    if (kvp.Value.ContainsKey(mod))
                    {
                        myLine += kvp.Value[mod] + "\t";
                    }
                    else
                    {
                        myLine += "0\t";
                    }
                }
                myOut.Add(myLine);
            }

            File.WriteAllLines(@"C:\Users\Michael Shortreed\Downloads\4-26-24 AD Brain DHAA Analysis\tableAllPeptidesMarkovich.txt", myOut);

        }

        [Test]
        public static void Junk6()
        {
            Dictionary<(int,string,string),List<string>> modPeptidesWithPositionLines = new Dictionary<(int, string, string), List<string>>();
            Dictionary<(int, string, string), List<string>> modPeptidesWithPositionFullSequences = new Dictionary<(int, string, string), List<string>>();
            Dictionary<(int, string, string), List<string>> modPeptidesWithPositionBaseSequences = new Dictionary<(int, string, string), List<string>>();
            Dictionary<string, List<string>> basePeptideToFullSequences = new Dictionary<string, List<string>>();

            using (StreamReader sr =
                   new StreamReader(
                       @"C:\Users\Michael Shortreed\Downloads\4-26-24 AD Brain DHAA Analysis\modPeptides180WithPosition.txt"))
            {
                bool continueReading = true;
                string line;
                bool firstLine = true;

                while ((line = sr.ReadLine()) != null && continueReading)
                {
                    if (!firstLine)
                    {
                        string[] fields = line.Split('\t');
                        string myKey = fields[0] + "\t" + fields[1] + "\t" + fields[2];

                        string baseSequence = fields[3];
                        while (baseSequence.Contains("[") && baseSequence.Contains("]"))
                        {
                            int firstOpenBracket = baseSequence.IndexOf('[');
                            int firstCloseBracket = baseSequence.IndexOf(']', firstOpenBracket);
                            if (firstCloseBracket != -1)
                            {
                                baseSequence = baseSequence.Remove(firstOpenBracket, firstCloseBracket - firstOpenBracket + 1);
                            }
                            else
                            {
                                break;
                            }
                        }
                        



                        if (modPeptidesWithPositionLines.ContainsKey((int.Parse(fields[0]), fields[1], fields[2])))
                        {
                            modPeptidesWithPositionLines[(int.Parse(fields[0]), fields[1], fields[2])].Add(line);
                        }
                        else
                        {
                            modPeptidesWithPositionLines.Add((int.Parse(fields[0]), fields[1], fields[2]), new List<string> { line });
                        }

                        if (modPeptidesWithPositionFullSequences.ContainsKey((int.Parse(fields[0]), fields[1],
                                fields[2])))
                        {
                            modPeptidesWithPositionFullSequences[(int.Parse(fields[0]), fields[1], fields[2])].Add(fields[3]);
                        }
                        else
                        {
                            modPeptidesWithPositionFullSequences.Add((int.Parse(fields[0]), fields[1], fields[2]), new List<string> { fields[3] });
                        }

                        if (modPeptidesWithPositionBaseSequences.ContainsKey((int.Parse(fields[0]), fields[1],
                                fields[2])))
                        {
                            modPeptidesWithPositionBaseSequences[(int.Parse(fields[0]), fields[1], fields[2])].Add(baseSequence);
                        }
                        else
                        {
                            modPeptidesWithPositionBaseSequences.Add((int.Parse(fields[0]), fields[1], fields[2]), new List<string> { baseSequence });
                        }

                        if (basePeptideToFullSequences.ContainsKey(baseSequence))
                        {
                            basePeptideToFullSequences[baseSequence].Add(fields[3]);
                        }
                        else
                        {
                            basePeptideToFullSequences.Add(baseSequence, new List<string> { fields[3] });
                        }
                    }
                    firstLine = false;
                }
            }

            foreach (var kvp in modPeptidesWithPositionBaseSequences)
            {
                modPeptidesWithPositionBaseSequences[kvp.Key] = kvp.Value.Distinct().ToList();
            }

            foreach (var kvp in basePeptideToFullSequences)
            {
                basePeptideToFullSequences[kvp.Key] = kvp.Value.Distinct().ToList();
            }
        }

        [Test]
        public static void Junk7()
        {
            Dictionary<(string, string), string> baseSequencefullSequencQvaluePepQvalueforPSMs = new Dictionary<(string, string), string>();

            using (StreamReader sr = new StreamReader(@"C:\Users\Michael Shortreed\Downloads\4-26-24 AD Brain DHAA Analysis\AllPSMsMarkovich.psmtsv"))
            {
                bool continueReading = true;
                string line;
                while ((line = sr.ReadLine()) != null && continueReading)
                {
                    string[] fields = line.Split('\t');
                    string baseSequence = fields[12];
                    string fullSequence = fields[13];
                    string targetDecoyContam = fields[38];
                    string qValue = fields[50];
                    string pepQValue = fields[55];
                    if (targetDecoyContam == "T" && double.TryParse(qValue, out double qValueDouble))
                    {
                        if (qValueDouble < 0.01 && !baseSequencefullSequencQvaluePepQvalueforPSMs.ContainsKey((baseSequence, fullSequence)))
                        {
                            baseSequencefullSequencQvaluePepQvalueforPSMs.Add((baseSequence, fullSequence), qValue + "\t" + pepQValue);
                        }

                    }

                    if (double.TryParse(qValue, out double qValueDouble2))
                    {
                        if (qValueDouble2 > 0.01)
                        {
                            continueReading = false;
                        }
                    }
                }
            }

            string psmFilePath =
                @"C:\Users\Michael Shortreed\Downloads\4-26-24 AD Brain DHAA Analysis\AllPeptidesMarkovich.psmtsv";
            List<PsmFromTsv> parsedPeptides = SpectrumMatchTsvReader.ReadPsmTsv(psmFilePath, out var warnings);
            parsedPeptides = parsedPeptides.Where(p => p.DecoyContamTarget.Contains("T")).ToList();

            List<string> interestingMods = new List<string> { "S[Common Biological:Phosphorylation on S]", "T[Common Biological:Phosphorylation on T]", "S[Less Common:Dehydroalanine on S]", "C[Less Common:Dehydroalanine on C]",
                "S[Custom:Homocys on S]", "C[Custom:Homocys on C]", "T[Custom:Homocys on T]", "T[Less Common:Dehydrobutyrine on T]", "S[Custom:DTT on S]", "C[Custom:DTT on C]", "C[Custom:DTT on T]", "T[Custom:Glutathione on T]",
                "S[Custom:Glutathione on S]", "C[Custom:Glutathione on C]", "S[Custom:TCEP on S]", "T[Custom:TCEP on T]", "C[Custom:TCEP on C]"  };

            Dictionary<(string, int), List<PsmFromTsv>> genePositionPeptides = new Dictionary<(string, int), List<PsmFromTsv>>();
            Dictionary<(string, int), List<string>> genePositionMod = new Dictionary<(string, int), List<string>>();
            Dictionary<string, List<PsmFromTsv>> noInterestingModsPeptideValues = new Dictionary<string, List<PsmFromTsv>>();

            foreach (var peptide in parsedPeptides)
            {
                if (baseSequencefullSequencQvaluePepQvalueforPSMs.ContainsKey((peptide.BaseSeq, peptide.FullSequence)))
                {
                    List<string> foundMods = interestingMods.Where(s => peptide.FullSequence.Contains(s)).ToList();
                    if (foundMods.Any())
                    {
                        foreach (var mod in foundMods)
                        {
                            string sequence = peptide.FullSequence;
                            string firstCharacter = mod.Substring(0, 1);
                            if (firstCharacter == "Y" || firstCharacter == "S" || firstCharacter == "T" ||
                                firstCharacter == "C")
                            {
                                sequence = sequence.Replace(mod, firstCharacter.ToLowerInvariant());
                            }

                            sequence = sequence.Replace("[I]", "");
                            sequence = sequence.Replace("[II]", "");
                            sequence = sequence.Replace("[III]", "");

                            //eliminate the remaning mods
                            while (sequence.Contains("[") && sequence.Contains("]"))
                            {
                                int firstOpenBracket = sequence.IndexOf('[');
                                int firstCloseBracket = sequence.IndexOf(']', firstOpenBracket);
                                if (firstCloseBracket != -1)
                                {
                                    sequence = sequence.Remove(firstOpenBracket, firstCloseBracket - firstOpenBracket + 1);
                                }
                                else
                                {
                                    break;
                                }
                            }

                            List<int> lowercaseLetterPositions = new List<int>();

                            for (int i = 0; i < sequence.Length; i++)
                            {
                                if (char.IsLower(sequence[i]))
                                {
                                    lowercaseLetterPositions.Add(i);
                                }
                            }

                            string firstAndLastAminoAcidPositionInProtein =
                                peptide.StartAndEndResiduesInProtein.Split('|')[0];
                            firstAndLastAminoAcidPositionInProtein =
                                firstAndLastAminoAcidPositionInProtein.Substring(1,
                                    firstAndLastAminoAcidPositionInProtein.Length - 2);
                            firstAndLastAminoAcidPositionInProtein =
                                firstAndLastAminoAcidPositionInProtein.Replace(" to ", "\t");
                            int[] startEnd = firstAndLastAminoAcidPositionInProtein.Split('\t').Select(int.Parse)
                                .ToArray();
                            lowercaseLetterPositions = lowercaseLetterPositions.Select(s => s + startEnd[0]).ToList();

                            string allGenesInPsm = peptide.GeneName;
                            string firstGene = "";
                            if (allGenesInPsm.Contains("|"))
                            {
                                string[] genes = allGenesInPsm.Split('|');

                                if (genes[0].Contains(":"))
                                {
                                    firstGene = genes[0].Split(':')[1];

                                    foreach (int position in lowercaseLetterPositions)
                                    {
                                        if (genePositionPeptides.ContainsKey((firstGene, position)))
                                        {
                                            genePositionPeptides[(firstGene, position)].Add(peptide);
                                            genePositionMod[(firstGene, position)]
                                                .Add(mod + "\t" + peptide.QValue + "\t" + peptide.PEP_QValue);

                                        }
                                        else
                                        {
                                            genePositionPeptides.Add((firstGene, position), new List<PsmFromTsv> { peptide });
                                            genePositionMod.Add((firstGene, position),
                                                new List<string> { mod + "\t" + peptide.QValue + "\t" + peptide.PEP_QValue });
                                        }
                                    }
                                }
                                else
                                {
                                    firstGene = genes[0];
                                    if (genes[0].Contains(":"))
                                    {
                                        firstGene = genes[0].Split(':')[1];
                                    }
                                    foreach (int position in lowercaseLetterPositions)
                                    {
                                        if (genePositionPeptides.ContainsKey((firstGene, position)))
                                        {
                                            genePositionPeptides[(firstGene, position)].Add(peptide);
                                            genePositionMod[(firstGene, position)]
                                                .Add(mod + "\t" + peptide.QValue + "\t" + peptide.PEP_QValue);

                                        }
                                        else
                                        {
                                            genePositionPeptides.Add((firstGene, position), new List<PsmFromTsv> { peptide });
                                            genePositionMod.Add((firstGene, position),
                                                new List<string> { mod + "\t" + peptide.QValue + "\t" + peptide.PEP_QValue });
                                        }
                                    }
                                }
                            }
                            else
                            {
                                firstGene = peptide.GeneName;
                                if (peptide.GeneName.Contains(":"))
                                {
                                    firstGene = peptide.GeneName.Split(':')[1];
                                }

                                foreach (int position in lowercaseLetterPositions)
                                {
                                    if (genePositionPeptides.ContainsKey((firstGene, position)))
                                    {
                                        genePositionPeptides[(firstGene, position)].Add(peptide);
                                        genePositionMod[(firstGene, position)]
                                            .Add(mod + "\t" + peptide.QValue + "\t" + peptide.PEP_QValue);

                                    }
                                    else
                                    {
                                        genePositionPeptides.Add((firstGene, position), new List<PsmFromTsv> { peptide });
                                        genePositionMod.Add((firstGene, position),
                                            new List<string> { mod + "\t" + peptide.QValue + "\t" + peptide.PEP_QValue });
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        if (noInterestingModsPeptideValues.ContainsKey(peptide.BaseSeq))
                        {
                            noInterestingModsPeptideValues[peptide.BaseSeq].Add(peptide);
                        }
                        else
                        {
                            noInterestingModsPeptideValues.Add(peptide.BaseSeq, new List<PsmFromTsv> { peptide });
                        }

                    }
                }
            }

            List<(string,string)> outList = new List<(string, string)>();
            List<string> myOut = new List<string>();
            myOut.Add("position" + "\t" + "protein accession" + "\t" + "gene" + "\t" + "base sequence" + "\t" + "full sequence" + "\t" + "modification" + "\t" + "Peptide Q-value" + "\t" + "Peptide PEP Q-Value" + "\t" + "PSM Q-value" + "\t" + "PSM PEP Q-Value");
            
            foreach (var kvp in genePositionPeptides)
            {
                if (kvp.Value.Count > 1)
                {
                    foreach (var psm in kvp.Value)
                    {
                        int index = kvp.Value.IndexOf(psm);
                        if (baseSequencefullSequencQvaluePepQvalueforPSMs.ContainsKey((psm.BaseSeq, psm.FullSequence)))
                        {
                            myOut.Add(kvp.Key.Item2 + "\t" + psm.ProteinAccession + "\t" + kvp.Key.Item1 + "\t" + psm.BaseSeq + "\t" + psm.FullSequence + "\t" + genePositionMod[kvp.Key][index] + "\t" + baseSequencefullSequencQvaluePepQvalueforPSMs[(psm.BaseSeq, psm.FullSequence)]);
                            outList.Add((psm.FullSequence, kvp.Key.Item2 + "\t" + psm.ProteinAccession + "\t" + kvp.Key.Item1 + "\t" + psm.BaseSeq + "\t" + psm.FullSequence + "\t" + genePositionMod[kvp.Key][index] + "\t" + baseSequencefullSequencQvaluePepQvalueforPSMs[(psm.BaseSeq, psm.FullSequence)]));
                            if (noInterestingModsPeptideValues.ContainsKey(psm.BaseSeq))
                            {
                                foreach (var psm2 in noInterestingModsPeptideValues[psm.BaseSeq])
                                {
                                    myOut.Add(kvp.Key.Item2 + "\t" + psm2.ProteinAccession + "\t" + kvp.Key.Item1 + "\t" + psm2.BaseSeq + "\t" + psm2.FullSequence + "\t" + "No Modifications of Interest" + "\t" + psm2.QValue + "\t" + psm2.PEP_QValue + "\t" + baseSequencefullSequencQvaluePepQvalueforPSMs[(psm2.BaseSeq, psm2.FullSequence)]);
                                    outList.Add((psm2.FullSequence, kvp.Key.Item2 + "\t" + psm2.ProteinAccession + "\t" + kvp.Key.Item1 + "\t" + psm2.BaseSeq + "\t" + psm2.FullSequence + "\t" + "No Modifications of Interest" + "\t" + psm2.QValue + "\t" + psm2.PEP_QValue + "\t" + baseSequencefullSequencQvaluePepQvalueforPSMs[(psm2.BaseSeq, psm2.FullSequence)]));
                                }
                                noInterestingModsPeptideValues.Remove(psm.BaseSeq);
                            }
                        }
                        else
                        {
                            //myOut.Add(kvp.Key.Item2 + "\t" + psm.ProteinAccession + "\t" + kvp.Key.Item1 + "\t" + psm.BaseSeq + "\t" + psm.FullSequence + "\t" + genePositionMod[kvp.Key][index]);
                            
                            //if (noInterestingModsPeptideValues.ContainsKey(psm.BaseSeq))
                            //{
                            //    foreach (var psm2 in noInterestingModsPeptideValues[psm.BaseSeq])
                            //    {
                            //        string psmQvalues = baseSequencefullSequencQvaluePepQvalueforPSMs[(psm2.BaseSeq, psm2.FullSequence)];
                            //        myOut.Add(kvp.Key.Item2 + "\t" + psm2.ProteinAccession + "\t" + kvp.Key.Item1 + "\t" + psm2.BaseSeq + "\t" + psm2.FullSequence + "\t" + "No Modifications of Interest" + "\t" + psm2.QValue + "\t" + psm2.PEP_QValue + "\t" + psmQvalues);

                            //    }
                            //    noInterestingModsPeptideValues.Remove(psm.BaseSeq);
                            //}
                        }

                    }
                }
            }

            Dictionary<string,string> quantDict = new Dictionary<string, string>();
            using (StreamReader sr =
                   new StreamReader(
                       @"C:\Users\Michael Shortreed\Downloads\4-26-24 AD Brain DHAA Analysis\QuantifiedPeptidesNormalizedMarkovich.tsv"))
            {
                string line;
                while ((line = sr.ReadLine()) != null)
                {
                    string[] fields = line.Split('\t');
                    string fullSequence = fields[0];
                    if (!quantDict.ContainsKey(fullSequence))
                    {
                        quantDict.Add(fullSequence, line);
                    }
                }
            }

            List<string> newOutList = new List<string>();
            newOutList.Add(myOut[0] + "\t" + quantDict["Sequence"]);
            foreach (var kvp in outList)
            {
                if (quantDict.ContainsKey(kvp.Item1))
                {
                    newOutList.Add(kvp.Item2 + "\t" + quantDict[kvp.Item1]);
                }
                else
                {
                    newOutList.Add(kvp.Item2);
                }
            }

            File.WriteAllLines(@"C:\Users\Michael Shortreed\Downloads\4-26-24 AD Brain DHAA Analysis\modPeptidesMarkovichWithPositionAndQuant_new.txt", newOutList);
        }


        [Test]
        public void SpectrumCount()
        {
            Assert.AreEqual(15, _mzSpectrumA.Size);
        }

        [Test]
        public void SpectrumFirstMZ()
        {
            Assert.AreEqual(328.73795, _mzSpectrumA.FirstX);
        }

        [Test]
        public void SpectrumLastMZ()
        {
            Assert.AreEqual(723.35345, _mzSpectrumA.LastX);
        }

        [Test]
        public void SpectrumBasePeakIntensity()
        {
            double basePeakIntensity = _mzSpectrumA.YofPeakWithHighestY.Value;

            Assert.AreEqual(122781408.0, basePeakIntensity);
        }

        [Test]
        public void SpectrumTIC()
        {
            double tic = _mzSpectrumA.SumOfAllY;

            Assert.AreEqual(843998894.0, tic);
        }

        [Test]
        public void SpectrumGetIntensityFirst()
        {
            Assert.AreEqual(81007096.0, _mzSpectrumA.YArray[0]);
        }

        [Test]
        public void SpectrumGetIntensityRandom()
        {
            Assert.AreEqual(44238040.0, _mzSpectrumA.YArray[6]);
        }

        [Test]
        public void SpectrumGetMassFirst()
        {
            Assert.AreEqual(328.73795, _mzSpectrumA.FirstX);
        }

        [Test]
        public void SpectrumGetMassRandom()
        {
            Assert.AreEqual(482.90393, _mzSpectrumA.XArray[6]);
        }

        [Test]
        public void SpectrumContainsPeak()
        {
            Assert.IsTrue(_mzSpectrumA.Size > 0);
        }

        [Test]
        public void SpectrumContainsPeakInRange()
        {
            Assert.AreEqual(1, _mzSpectrumA.NumPeaksWithinRange(448.23987 - 0.001, 448.23987 + 0.001));
        }

        // Tests in current implementation
        [Test]
        public void SpectrumContainsPeakInRangeEnd()
        {
            Assert.AreEqual(1, _mzSpectrumA.NumPeaksWithinRange(448.23987 - 0.001, 448.23987));
        }

        [Test]
        public void SpectrumContainsPeakInRangeStart()
        {
            Assert.AreEqual(1, _mzSpectrumA.NumPeaksWithinRange(448.23987, 448.23987 + 0.001));
        }

        [Test]
        public void SpectrumContainsPeakInRangeStartEnd()
        {
            Assert.AreEqual(1, _mzSpectrumA.NumPeaksWithinRange(448.23987, 448.23987));
        }

        [Test]
        public void SpectrumDoesntContainPeakInRange()
        {
            Assert.AreEqual(0, _mzSpectrumA.NumPeaksWithinRange(603.4243 - 0.001, 603.4243 + 0.001));
        }

        [Test]
        public void SpectrumMassRange()
        {
            MzRange range = new MzRange(328.73795, 723.35345);

            Assert.AreEqual(0, _mzSpectrumA.Range.Minimum - range.Minimum, 1e-9);
            Assert.AreEqual(0, _mzSpectrumA.Range.Maximum - range.Maximum, 1e-9);
        }

        [Test]
        public void SpectrumFilterCount()
        {
            var filteredMzSpectrum = _mzSpectrumA.FilterByY(28604417, 28604419);

            Assert.AreEqual(1, filteredMzSpectrum.Count());
        }

        [Test]
        public void FilterByNumberOfMostIntenseTest()
        {
            Assert.AreEqual(5, _mzSpectrumA.FilterByNumberOfMostIntense(5).Count());
        }

        [Test]
        public void FilterByNumberOfMostIntenseRobTest()
        {
            double[] x = new double[] { 50, 60, 70, 147.0764, 257.1244, 258.127, 275.135 };
            double[] y = new double[] { 1, 1, 1, 1, 1, 1, 1 };
            MzSpectrum spectrum = new MzSpectrum(x, y, false);
            Assert.AreEqual(7, spectrum.FilterByNumberOfMostIntense(200).Count());
        }

        [Test]
        public void GetBasePeak()
        {
            Assert.AreEqual(122781408.0, _mzSpectrumA.YofPeakWithHighestY);
        }

        [Test]
        public void GetClosestPeak()
        {
            Assert.AreEqual(448.23987, _mzSpectrumA.GetClosestPeakXvalue(448));
            Assert.AreEqual(447.73849, _mzSpectrumA.GetClosestPeakXvalue(447.9));
        }

        [Test]
        public void Extract()
        {
            Assert.AreEqual(3, _mzSpectrumA.Extract(500, 600).Count());
        }

        [Test]
        public void CorrectOrder()
        {
            _mzSpectrumA = new MzSpectrum(new double[] { 5, 6, 7 }, new double[] { 1, 2, 3 }, false);
            Assert.IsTrue(_mzSpectrumA.FilterByNumberOfMostIntense(2).First().Mz < _mzSpectrumA.FilterByNumberOfMostIntense(2).ToList()[1].Mz);
        }

        [Test]
        public void TestFunctionToX()
        {
            _mzSpectrumA.ReplaceXbyApplyingFunction(b => -1);
            Assert.AreEqual(-1, _mzSpectrumA.XArray[0]);
        }

        [Test]
        public void TestGetClosestPeakXValue()
        {
            Assert.AreEqual(447.73849, _mzSpectrumA.GetClosestPeakXvalue(447.73849));
            Assert.AreEqual(447.73849, _mzSpectrumA.GetClosestPeakXvalue(447));
            Assert.IsNull(new MzSpectrum(new double[0], new double[0], false).GetClosestPeakXvalue(1));
        }

        [Test]
        public void TestDotProduct()
        {
            double[] array1 = { 1 };
            double[] array2 = { 2 };
            double[] array3 = { 1, 2 };
            double[] array4 = { 1, 1 };

            MzSpectrum spec1 = new MzSpectrum(array1, array1, false);
            MzSpectrum spec2 = new MzSpectrum(array2, array1, false);
            MzSpectrum spec3 = new MzSpectrum(array3, array4, false);
            Tolerance tolerance = new PpmTolerance(10);

            Assert.AreEqual(spec1.CalculateDotProductSimilarity(spec3, tolerance), spec3.CalculateDotProductSimilarity(spec1, tolerance)); //comparison side shouldn't matter
            Assert.AreEqual(spec1.CalculateDotProductSimilarity(spec2, tolerance), 0); //orthogonal spectra give a score of zero
            Assert.AreEqual(spec2.CalculateDotProductSimilarity(spec2, tolerance), 1); //identical spectra give a score of 1
            Assert.IsTrue(tolerance.Within(spec3.CalculateDotProductSimilarity(spec2, tolerance), Math.Cos(Math.PI / 4)));
        }

        [Test]
        public void TestNumPeaksWithinRange()
        {
            double[] xArray = { 1, 2, 3, 4, 5, 6, 7 };
            double[] yArray = { 1, 2, 1, 5, 1, 2, 1 };
            var thisSpectrum = new MzSpectrum(xArray, yArray, false);

            Assert.AreEqual(7, thisSpectrum.NumPeaksWithinRange(double.MinValue, double.MaxValue));

            Assert.AreEqual(7, thisSpectrum.NumPeaksWithinRange(1, 7));

            Assert.AreEqual(1, thisSpectrum.NumPeaksWithinRange(1, 1));

            Assert.AreEqual(2, thisSpectrum.NumPeaksWithinRange(1, 2));

            Assert.AreEqual(2, thisSpectrum.NumPeaksWithinRange(0.001, 2.999));

            Assert.AreEqual(1, thisSpectrum.NumPeaksWithinRange(0, 1.5));

            Assert.AreEqual(1, thisSpectrum.NumPeaksWithinRange(6.5, 8));

            Assert.AreEqual(3, thisSpectrum.NumPeaksWithinRange(3, 5));

            Assert.AreEqual(2, thisSpectrum.NumPeaksWithinRange(3.5, 5.5));

            Assert.AreEqual(1, thisSpectrum.NumPeaksWithinRange(7, 8));

            Assert.AreEqual(0, thisSpectrum.NumPeaksWithinRange(8, 9));

            Assert.AreEqual(0, thisSpectrum.NumPeaksWithinRange(-2, -1));

            Assert.AreEqual("[1 to 7] m/z (Peaks 7)", thisSpectrum.ToString());

            //Assert.AreEqual(7, thisSpectrum.FilterByNumberOfMostIntense(7).Size);
            //Assert.AreEqual(1, thisSpectrum.FilterByNumberOfMostIntense(1).Size);
            //Assert.AreEqual(4, thisSpectrum.FilterByNumberOfMostIntense(1).FirstX);

            //Assert.AreEqual(2, thisSpectrum.FilterByNumberOfMostIntense(3).FirstX);

            //Assert.AreEqual(0, thisSpectrum.FilterByNumberOfMostIntense(0).Size);

            //Assert.AreEqual(2, thisSpectrum.WithRangeRemoved(2, 6).Size);
            //Assert.AreEqual(0, thisSpectrum.WithRangeRemoved(0, 100).Size);

            //Assert.AreEqual(6, thisSpectrum.WithRangeRemoved(7, 100).Size);

            //Assert.AreEqual(1, thisSpectrum.WithRangeRemoved(new DoubleRange(double.MinValue, 6)).Size);

            List<DoubleRange> xRanges = new List<DoubleRange>
            {
                new DoubleRange(2, 5),
                new DoubleRange(3, 6)
            };
            //Assert.AreEqual(2, thisSpectrum.WithRangesRemoved(xRanges).Size);

            //Assert.AreEqual(3, thisSpectrum.Extract(new DoubleRange(4.5, 10)).Size);

            //Assert.AreEqual(2, thisSpectrum.FilterByY(new DoubleRange(1.5, 2.5)).Size);

            //Assert.AreEqual(3, thisSpectrum.FilterByY(1.5, double.MaxValue).Size);

            //Assert.AreEqual(2, thisSpectrum.ApplyFunctionToX(b => b * 2).FirstX);

            Assert.AreEqual(1, thisSpectrum.GetClosestPeakXvalue(-100));

            Assert.AreEqual(7, thisSpectrum.GetClosestPeakXvalue(6.6));

            Assert.AreEqual(7, thisSpectrum.GetClosestPeakXvalue(7));

            Assert.AreEqual(7, thisSpectrum.GetClosestPeakXvalue(8));
        }

        [Test]
        public void TestEqualsAndHashCode()
        {
            // identical spectra, x and y arrays deep copied
            MzSpectrum identicalSpectrum = new(_mzSpectrumA.XArray, _mzSpectrumA.YArray, true);
            Assert.AreEqual(identicalSpectrum.GetHashCode(), _mzSpectrumA.GetHashCode());
            Assert.IsTrue(identicalSpectrum.Equals(_mzSpectrumA));
            Assert.IsTrue(identicalSpectrum.Equals((object)_mzSpectrumA));

            // changed x value
            identicalSpectrum.XArray[1] += 10;
            Assert.IsFalse(identicalSpectrum.Equals(_mzSpectrumA));
            Assert.IsFalse(identicalSpectrum.Equals((object)_mzSpectrumA));
            Assert.AreNotEqual(identicalSpectrum.GetHashCode(), _mzSpectrumA.GetHashCode());
            identicalSpectrum.XArray[1] -= 10;

            // changed y value
            identicalSpectrum.YArray[1] += 10;
            Assert.IsFalse(identicalSpectrum.Equals(_mzSpectrumA));
            Assert.IsFalse(identicalSpectrum.Equals((object)_mzSpectrumA));
            Assert.AreNotEqual(identicalSpectrum.GetHashCode(), _mzSpectrumA.GetHashCode());

            Assert.That(!_mzSpectrumA.Equals(null));
            Assert.That(!_mzSpectrumA.Equals((object)null));
            Assert.That(!_mzSpectrumA.Equals(2));
            Assert.That(!_mzSpectrumA.Equals((object)2));
        }
    }
}