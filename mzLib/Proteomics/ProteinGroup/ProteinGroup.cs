using Easy.Common.Extensions;
using MassSpectrometry;
using MzLibUtil;
using Omics;
using Omics.Modifications;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using Omics.BioPolymerGroup;

namespace Proteomics.ProteinGroup
{
    public class ProteinGroup : IBioPolymerGroup
    {
        /// <summary>
        /// Maximum length for string fields in output. Strings longer than this will be truncated.
        /// Set to 0 or negative to disable truncation.
        /// </summary>
        public static int MaxStringLength { get; set; } = 32000;

        public ProteinGroup(HashSet<IBioPolymer> proteins, HashSet<IBioPolymerWithSetMods> peptides,
            HashSet<IBioPolymerWithSetMods> uniquePeptides)
        {
            BioPolymers = proteins;
            ListOfBioPolymersOrderedByAccession = BioPolymers.OrderBy(p => p.Accession).ToList();
            BioPolymerGroupName = string.Join("|", ListOfBioPolymersOrderedByAccession.Select(p => p.Accession));
            AllBioPolymersWithSetMods = peptides;
            UniqueBioPolymersWithSetMods = uniquePeptides;
            AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch>();
            SequenceCoverageFraction = new List<double>();
            SequenceCoverageDisplayList = new List<string>();
            SequenceCoverageDisplayListWithMods = new List<string>();
            FragmentSequenceCoverageDisplayList = new List<string>();
            BioPolymerGroupScore = 0;
            BestBioPolymerWithSetModsScore = 0;
            QValue = 0;
            IsDecoy = false;
            IsContaminant = false;
            ModsInfo = new List<string>();

            // if any of the proteins in the protein group are decoys, the protein group is a decoy
            foreach (var protein in proteins)
            {
                if (protein.IsDecoy)
                {
                    IsDecoy = true;
                    break;
                }

                if (protein.IsContaminant)
                {
                    IsContaminant = true;
                    break;
                }
            }
        }

        #region IBioPolymerGroup Implementation

        public bool IsDecoy { get; }

        public bool IsContaminant { get; }

        public List<ISampleInfo> SamplesForQuantification { get; set; }

        public HashSet<IBioPolymer> BioPolymers { get; set; }

        public string BioPolymerGroupName { get; private set; }

        public double BioPolymerGroupScore { get; set; }

        public HashSet<IBioPolymerWithSetMods> AllBioPolymersWithSetMods { get; set; }

        public HashSet<IBioPolymerWithSetMods> UniqueBioPolymersWithSetMods { get; set; }

        public HashSet<ISpectralMatch> AllPsmsBelowOnePercentFDR { get; set; }

        public double QValue { get; set; }

        public double BestBioPolymerWithSetModsQValue { get; set; }

        public double BestBioPolymerWithSetModsScore { get; set; }

        public List<string> ModsInfo { get; private set; }

        public Dictionary<ISampleInfo, double> IntensitiesBySample { get; set; }

        public List<IBioPolymer> ListOfBioPolymersOrderedByAccession { get; private set; }

        #endregion

        #region Legacy Properties (for backward compatibility)

        /// <summary>
        /// Legacy property. Use <see cref="SamplesForQuantification"/> instead.
        /// Returns only SpectraFileInfo samples for label-free compatibility.
        /// </summary>
        public List<SpectraFileInfo> FilesForQuantification
        {
            get => SamplesForQuantification?.OfType<SpectraFileInfo>().ToList();
            set => SamplesForQuantification = value?.Cast<ISampleInfo>().ToList();
        }

        /// <summary>
        /// Legacy property. Use <see cref="BioPolymers"/> instead.
        /// </summary>
        public HashSet<IBioPolymer> Proteins
        {
            get => BioPolymers;
            set => BioPolymers = value;
        }

        /// <summary>
        /// Legacy property. Use <see cref="BioPolymerGroupName"/> instead.
        /// </summary>
        public string ProteinGroupName => BioPolymerGroupName;

        /// <summary>
        /// Legacy property. Use <see cref="BioPolymerGroupScore"/> instead.
        /// </summary>
        public double ProteinGroupScore
        {
            get => BioPolymerGroupScore;
            set => BioPolymerGroupScore = value;
        }

        /// <summary>
        /// Legacy property. Use <see cref="AllBioPolymersWithSetMods"/> instead.
        /// </summary>
        public HashSet<IBioPolymerWithSetMods> AllPeptides
        {
            get => AllBioPolymersWithSetMods;
            set => AllBioPolymersWithSetMods = value;
        }

        /// <summary>
        /// Legacy property. Use <see cref="UniqueBioPolymersWithSetMods"/> instead.
        /// </summary>
        public HashSet<IBioPolymerWithSetMods> UniquePeptides
        {
            get => UniqueBioPolymersWithSetMods;
            set => UniqueBioPolymersWithSetMods = value;
        }

        /// <summary>
        /// Legacy property. Use <see cref="BestBioPolymerWithSetModsQValue"/> instead.
        /// </summary>
        public double BestPeptideQValue
        {
            get => BestBioPolymerWithSetModsQValue;
            set => BestBioPolymerWithSetModsQValue = value;
        }

        /// <summary>
        /// Legacy property. Use <see cref="BestBioPolymerWithSetModsScore"/> instead.
        /// </summary>
        public double BestPeptideScore
        {
            get => BestBioPolymerWithSetModsScore;
            set => BestBioPolymerWithSetModsScore = value;
        }

        /// <summary>
        /// Legacy property. Use <see cref="IntensitiesBySample"/> instead.
        /// Returns only SpectraFileInfo entries for label-free compatibility.
        /// </summary>
        public Dictionary<SpectraFileInfo, double> IntensitiesByFile
        {
            get => IntensitiesBySample?
                .Where(kvp => kvp.Key is SpectraFileInfo)
                .ToDictionary(kvp => (SpectraFileInfo)kvp.Key, kvp => kvp.Value);
            set => IntensitiesBySample = value?.ToDictionary(kvp => (ISampleInfo)kvp.Key, kvp => kvp.Value);
        }

        #endregion

        #region Additional Properties

        public List<double> SequenceCoverageFraction { get; private set; }

        public List<string> SequenceCoverageDisplayList { get; private set; }

        public List<string> SequenceCoverageDisplayListWithMods { get; private set; }

        public List<string> FragmentSequenceCoverageDisplayList { get; private set; }

        public int CumulativeTarget { get; set; }

        public int CumulativeDecoy { get; set; }

        public bool DisplayModsOnPeptides { get; set; }

        private string UniquePeptidesOutput;
        private string SharedPeptidesOutput;

        #endregion

        #region Methods

        public void GetIdentifiedPeptidesOutput(List<SilacLabel> labels)
        {
            var SharedPeptides = AllBioPolymersWithSetMods.Except(UniqueBioPolymersWithSetMods);
            if (labels != null)
            {
                // SILAC handling - currently commented out
            }
            else
            {
                if (!DisplayModsOnPeptides)
                {
                    UniquePeptidesOutput =
                        TruncateString(string.Join("|",
                            UniqueBioPolymersWithSetMods.Select(p => p.BaseSequence).Distinct()));
                    SharedPeptidesOutput =
                        TruncateString(string.Join("|",
                            SharedPeptides.Select(p => p.BaseSequence).Distinct()));
                }
                else
                {
                    UniquePeptidesOutput =
                        TruncateString(string.Join("|",
                            UniqueBioPolymersWithSetMods.Select(p => p.FullSequence).Distinct()));
                    SharedPeptidesOutput =
                        TruncateString(string.Join("|",
                            SharedPeptides.Select(p => p.FullSequence).Distinct()));
                }
            }
        }

        public string GetTabSeparatedHeader()
        {
            var sb = new StringBuilder();
            sb.Append("Protein Accession" + '\t');
            sb.Append("Gene" + '\t');
            sb.Append("Organism" + '\t');
            sb.Append("Protein Full Name" + '\t');
            sb.Append("Protein Unmodified Mass" + '\t');
            sb.Append("Number of Proteins in Group" + '\t');
            sb.Append("Unique Peptides" + '\t');
            sb.Append("Shared Peptides" + '\t');
            sb.Append("Number of Peptides" + '\t');
            sb.Append("Number of Unique Peptides" + '\t');
            sb.Append("Sequence Coverage Fraction" + '\t');
            sb.Append("Sequence Coverage" + '\t');
            sb.Append("Sequence Coverage with Mods" + '\t');
            sb.Append("Fragment Sequence Coverage" + '\t');
            sb.Append("Modification Info List" + "\t");

            if (SamplesForQuantification != null)
            {
                // Check if we have label-free (SpectraFileInfo) or isobaric (IsobaricQuantSampleInfo) samples
                var spectraFiles = SamplesForQuantification.OfType<SpectraFileInfo>().ToList();
                var isobaricSamples = SamplesForQuantification.OfType<IsobaricQuantSampleInfo>().ToList();

                if (spectraFiles.Any())
                {
                    // Label-free header generation
                    bool unfractionated = spectraFiles.Select(p => p.Fraction).Distinct().Count() == 1;
                    bool conditionsUndefined = spectraFiles.All(p => string.IsNullOrEmpty(p.Condition));
                    bool silacExperimentalDesign = spectraFiles.Any(p => !File.Exists(p.FullFilePathWithExtension));

                    foreach (var sampleGroup in spectraFiles.GroupBy(p => p.Condition))
                    {
                        foreach (var sample in sampleGroup.GroupBy(p => p.BiologicalReplicate).OrderBy(p => p.Key))
                        {
                            if ((conditionsUndefined && unfractionated) || silacExperimentalDesign)
                            {
                                sb.Append("Intensity_" + sample.First().FilenameWithoutExtension + "\t");
                            }
                            else
                            {
                                sb.Append("Intensity_" + sample.First().Condition + "_" +
                                          (sample.First().BiologicalReplicate + 1) + "\t");
                            }
                        }
                    }
                }
                else if (isobaricSamples.Any())
                {
                    // Isobaric header generation - group by file, then by channel
                    foreach (var fileGroup in isobaricSamples.GroupBy(p => p.FullFilePathWithExtension).OrderBy(g => g.Key))
                    {
                        foreach (var sample in fileGroup.OrderBy(p => p.ChannelLabel))
                        {
                            sb.Append($"Intensity_{Path.GetFileNameWithoutExtension(sample.FullFilePathWithExtension)}_{sample.ChannelLabel}\t");
                        }
                    }
                }
            }

            sb.Append("Number of PSMs" + '\t');
            sb.Append("Protein Decoy/Contaminant/Target" + '\t');
            sb.Append("Protein Cumulative Target" + '\t');
            sb.Append("Protein Cumulative Decoy" + '\t');
            sb.Append("Protein QValue" + '\t');
            sb.Append("Best Peptide Score" + '\t');
            sb.Append("Best Peptide Notch QValue");
            return sb.ToString();
        }

        public override string ToString()
        {
            var sb = new StringBuilder();

            sb.Append(BioPolymerGroupName);
            sb.Append("\t");

            sb.Append(TruncateString(string.Join("|",
                ListOfBioPolymersOrderedByAccession.Select(p => p.GeneNames.Select(x => x.Item2).FirstOrDefault()))));
            sb.Append("\t");

            sb.Append(TruncateString(string.Join("|",
                ListOfBioPolymersOrderedByAccession.Select(p => p.Organism).Distinct())));
            sb.Append("\t");

            sb.Append(TruncateString(string.Join("|",
                ListOfBioPolymersOrderedByAccession.Select(p => p.FullName).Distinct())));
            sb.Append("\t");

            var sequences = ListOfBioPolymersOrderedByAccession.Select(p => p.BaseSequence).Distinct();
            List<double> masses = new List<double>();
            foreach (var sequence in sequences)
            {
                try
                {
                    masses.Add(new Proteomics.AminoAcidPolymer.Peptide(sequence).MonoisotopicMass);
                }
                catch (Exception)
                {
                    masses.Add(double.NaN);
                }
            }

            sb.Append(TruncateString(string.Join("|", masses)));
            sb.Append("\t");

            sb.Append("" + BioPolymers.Count);
            sb.Append("\t");

            if (UniquePeptidesOutput != null)
            {
                sb.Append(TruncateString(UniquePeptidesOutput));
            }
            sb.Append("\t");

            if (SharedPeptidesOutput != null)
            {
                sb.Append(TruncateString(SharedPeptidesOutput));
            }
            sb.Append("\t");

            if (!DisplayModsOnPeptides)
            {
                sb.Append("" + AllBioPolymersWithSetMods.Select(p => p.BaseSequence).Distinct().Count());
            }
            else
            {
                sb.Append("" + AllBioPolymersWithSetMods.Select(p => p.FullSequence).Distinct().Count());
            }
            sb.Append("\t");

            if (!DisplayModsOnPeptides)
            {
                sb.Append("" + UniqueBioPolymersWithSetMods.Select(p => p.BaseSequence).Distinct().Count());
            }
            else
            {
                sb.Append("" + UniqueBioPolymersWithSetMods.Select(p => p.FullSequence).Distinct().Count());
            }
            sb.Append("\t");

            sb.Append(TruncateString(string.Join("|",
                SequenceCoverageFraction.Select(p => string.Format("{0:0.#####}", p)))));
            sb.Append("\t");

            sb.Append(TruncateString(string.Join("|", SequenceCoverageDisplayList)));
            sb.Append("\t");

            sb.Append(TruncateString(string.Join("|", SequenceCoverageDisplayListWithMods)));
            sb.Append("\t");

            sb.Append(TruncateString(string.Join("|", FragmentSequenceCoverageDisplayList)));
            sb.Append("\t");

            sb.Append(TruncateString(string.Join("|", ModsInfo)));
            sb.Append("\t");

            // Output intensities
            if (IntensitiesBySample != null && SamplesForQuantification != null)
            {
                var spectraFiles = SamplesForQuantification.OfType<SpectraFileInfo>().ToList();
                var isobaricSamples = SamplesForQuantification.OfType<IsobaricQuantSampleInfo>().ToList();

                if (spectraFiles.Any())
                {
                    // Label-free intensity output
                    foreach (var sampleGroup in spectraFiles.GroupBy(p => p.Condition))
                    {
                        foreach (var sample in sampleGroup.GroupBy(p => p.BiologicalReplicate).OrderBy(p => p.Key))
                        {
                            double summedIntensity = sample.Sum(file =>
                                IntensitiesBySample.TryGetValue(file, out var intensity) ? intensity : 0);

                            if (summedIntensity > 0)
                            {
                                sb.Append(summedIntensity);
                            }
                            sb.Append("\t");
                        }
                    }
                }
                else if (isobaricSamples.Any())
                {
                    // Isobaric intensity output
                    foreach (var fileGroup in isobaricSamples.GroupBy(p => p.FullFilePathWithExtension).OrderBy(g => g.Key))
                    {
                        foreach (var sample in fileGroup.OrderBy(p => p.ChannelLabel))
                        {
                            if (IntensitiesBySample.TryGetValue(sample, out var intensity) && intensity > 0)
                            {
                                sb.Append(intensity);
                            }
                            sb.Append("\t");
                        }
                    }
                }
            }

            sb.Append("" + AllPsmsBelowOnePercentFDR.Count);
            sb.Append("\t");

            if (IsDecoy)
            {
                sb.Append("D");
            }
            else if (IsContaminant)
            {
                sb.Append("C");
            }
            else
            {
                sb.Append("T");
            }

            sb.Append("\t");
            sb.Append(CumulativeTarget);
            sb.Append("\t");
            sb.Append(CumulativeDecoy);
            sb.Append("\t");
            sb.Append(QValue);
            sb.Append("\t");
            sb.Append(BestBioPolymerWithSetModsScore);
            sb.Append("\t");
            sb.Append(BestBioPolymerWithSetModsQValue);

            return sb.ToString();
        }

        public void Score()
        {
            BioPolymerGroupScore = AllPsmsBelowOnePercentFDR.GroupBy(p => p.BaseSequence)
                .Select(p => p.Select(x => x.Score).Max()).Sum();
        }

        public void MergeWith(IBioPolymerGroup otherBioPolymerGroup)
        {
            if (otherBioPolymerGroup is ProteinGroup other)
            {
                MergeProteinGroupWith(other);
            }
        }

        public void MergeProteinGroupWith(ProteinGroup other)
        {
            this.BioPolymers.UnionWith(other.BioPolymers);
            this.AllBioPolymersWithSetMods.UnionWith(other.AllBioPolymersWithSetMods);
            this.UniqueBioPolymersWithSetMods.UnionWith(other.UniqueBioPolymersWithSetMods);
            this.AllPsmsBelowOnePercentFDR.UnionWith(other.AllPsmsBelowOnePercentFDR);
            other.BioPolymerGroupScore = 0;

            ListOfBioPolymersOrderedByAccession = BioPolymers.OrderBy(p => p.Accession).ToList();
            BioPolymerGroupName = string.Join("|", ListOfBioPolymersOrderedByAccession.Select(p => p.Accession));
        }

        public IBioPolymerGroup ConstructSubsetBioPolymerGroup(string fullFilePath, List<SilacLabel>? silacLabels = null)
        {
            return ConstructSubsetProteinGroup(fullFilePath, silacLabels);
        }

        public ProteinGroup ConstructSubsetProteinGroup(string fullFilePath, List<SilacLabel>? silacLabels = null)
        {
            var allPsmsForThisFile =
                new HashSet<ISpectralMatch>(
                    AllPsmsBelowOnePercentFDR.Where(p => p.FullFilePath.Equals(fullFilePath)));
            var allPeptidesForThisFile =
                new HashSet<IBioPolymerWithSetMods>(
                    allPsmsForThisFile.SelectMany(p => p.GetIdentifiedBioPolymersWithSetMods()));
            var allUniquePeptidesForThisFile =
                new HashSet<IBioPolymerWithSetMods>(UniqueBioPolymersWithSetMods.Intersect(allPeptidesForThisFile));

            ProteinGroup subsetPg = new ProteinGroup(BioPolymers, allPeptidesForThisFile, allUniquePeptidesForThisFile)
            {
                AllPsmsBelowOnePercentFDR = allPsmsForThisFile,
                DisplayModsOnPeptides = DisplayModsOnPeptides
            };

            if (SamplesForQuantification != null)
            {
                // Find matching sample(s) for this file path
                var matchingSamples = SamplesForQuantification
                    .Where(p => p.FullFilePathWithExtension == fullFilePath)
                    .ToList();

                subsetPg.SamplesForQuantification = matchingSamples;

                if (IntensitiesBySample != null)
                {
                    subsetPg.IntensitiesBySample = IntensitiesBySample
                        .Where(kvp => matchingSamples.Contains(kvp.Key))
                        .ToDictionary(kvp => kvp.Key, kvp => kvp.Value);
                }
            }

            return subsetPg;
        }

        #endregion

        #region Equality

        public bool Equals(IBioPolymerGroup? other)
        {
            if (other is null) return false;
            if (ReferenceEquals(this, other)) return true;
            return BioPolymerGroupName == other.BioPolymerGroupName;
        }

        public bool Equals(ProteinGroup? grp)
        {
            if (grp == null) return false;
            return BioPolymerGroupName == grp.BioPolymerGroupName;
        }

        public override bool Equals(object? obj)
        {
            if (obj is ProteinGroup pg) return Equals(pg);
            if (obj is IBioPolymerGroup bg) return Equals(bg);
            return false;
        }

        public override int GetHashCode()
        {
            return BioPolymerGroupName?.GetHashCode() ?? 0;
        }

        #endregion

        #region Private Helpers

        private static string TruncateString(string? input)
        {
            if (string.IsNullOrEmpty(input))
                return string.Empty;

            if (MaxStringLength <= 0 || input.Length <= MaxStringLength)
                return input;

            return input.Substring(0, MaxStringLength);
        }

        #endregion
    }
}