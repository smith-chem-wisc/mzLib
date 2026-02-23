using Easy.Common.Extensions;
using MassSpectrometry;
using Omics.Modifications;
using System.Text;

namespace Omics.BioPolymerGroup
{
    /// <summary>
    /// Represents a group of related biopolymers (e.g., proteins, RNA sequences) that share 
    /// identified peptide or oligonucleotide sequences. Groups are formed during protein/gene 
    /// inference when multiple biopolymers cannot be distinguished based on the identified sequences.
    /// 
    /// This class provides:
    /// <list type="bullet">
    ///   <item><description>Sequence coverage calculation at both peptide-level and fragment-level</description></item>
    ///   <item><description>Quantification support for label-free (spectral counting) and isobaric (TMT/iTRAQ) methods</description></item>
    ///   <item><description>Modification occupancy statistics</description></item>
    ///   <item><description>FDR calculation support via cumulative target/decoy counting</description></item>
    ///   <item><description>Tab-separated output formatting for results files</description></item>
    /// </list>
    /// 
    /// Fragment-level coverage is only calculated when PSMs implement <see cref="IHasSequenceCoverageFromFragments"/>.
    /// </summary>
    public class BioPolymerGroup : IBioPolymerGroup
    {
        /// <summary>
        /// Maximum length for string fields in output. Strings exceeding this length will be truncated.
        /// 
        /// Default is 32,000 characters, which is slightly below Excel's cell limit of 32,767 characters.
        /// This ensures output files can be opened in Excel without data truncation or corruption.
        /// Set to 0 or negative to disable truncation (useful for programmatic processing where 
        /// Excel compatibility is not required).
        /// </summary>
        /// <remarks>
        /// Excel specification: A cell can contain up to 32,767 characters.
        /// See: https://support.microsoft.com/en-us/office/excel-specifications-and-limits
        /// </remarks>
        public static int MaxStringLength { get; set; } = 32000;

        /// <summary>
        /// Creates a new biopolymer group from the specified biopolymers and identified sequences.
        /// </summary>
        /// <param name="bioPolymers">Set of biopolymers (e.g., proteins, RNA) that belong to this group.
        /// These are typically indistinguishable based on the identified sequences.</param>
        /// <param name="bioPolymersWithSetMods">All identified sequences with modifications for this group,
        /// including sequences shared with other groups.</param>
        /// <param name="uniqueBioPolymersWithSetMods">Sequences with modifications that are unique to this group
        /// and not shared with any other biopolymer group.</param>
        public BioPolymerGroup(HashSet<IBioPolymer> bioPolymers, HashSet<IBioPolymerWithSetMods> bioPolymersWithSetMods,
            HashSet<IBioPolymerWithSetMods> uniqueBioPolymersWithSetMods)
        {
            BioPolymers = bioPolymers;
            ListOfBioPolymersOrderedByAccession = BioPolymers.OrderBy(p => p.Accession).ToList();
            BioPolymerGroupName = string.Join("|", ListOfBioPolymersOrderedByAccession.Select(p => p.Accession));
            AllBioPolymersWithSetMods = bioPolymersWithSetMods;
            UniqueBioPolymersWithSetMods = uniqueBioPolymersWithSetMods;
            AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch>();
            BioPolymerGroupScore = 0;
            BestBioPolymerWithSetModsScore = 0;
            QValue = 0;
            IsDecoy = false;
            IsContaminant = false;

            // if any of the biopolymers in the group are decoys, the group is a decoy
            foreach (var bioPolymer in bioPolymers)
            {
                if (bioPolymer.IsDecoy)
                {
                    IsDecoy = true;
                }

                if (bioPolymer.IsContaminant)
                {
                    IsContaminant = true;
                }

                // If both are true, we can break early
                if (IsDecoy && IsContaminant)
                {
                    break;
                }
            }
        }

        #region IBioPolymerGroup Implementation

        /// <summary>
        /// True if this group contains any decoy biopolymers, used for FDR estimation.
        /// </summary>
        public bool IsDecoy { get; }

        /// <summary>
        /// True if this group contains any biopolymers marked as contaminants.
        /// </summary>
        public bool IsContaminant { get; }

        /// <summary>
        /// List of samples that contribute quantification data for this group.
        /// Supports both <see cref="SpectraFileInfo"/> (label-free) and <see cref="IsobaricQuantSampleInfo"/> (TMT/iTRAQ).
        /// </summary>
        public List<ISampleInfo> SamplesForQuantification { get; set; }

        /// <summary>
        /// Set of all biopolymers (e.g., proteins, RNA sequences) that belong to this group.
        /// These biopolymers are indistinguishable based on the identified sequences.
        /// </summary>
        public HashSet<IBioPolymer> BioPolymers { get; set; }

        /// <summary>
        /// Display name for the biopolymer group, derived from the pipe-delimited accessions of member biopolymers.
        /// Used as the primary identity key for equality comparisons via <see cref="IEquatable{T}"/>.
        /// </summary>
        public string BioPolymerGroupName { get; private set; }

        /// <summary>
        /// Aggregated confidence score for the group. Higher values indicate higher confidence.
        /// 
        /// Computed by <see cref="BioPolymerGroupExtensions.Score"/> as the sum of the best (highest) 
        /// score for each unique base sequence among the PSMs in <see cref="AllPsmsBelowOnePercentFDR"/>.
        /// This ensures each unique peptide/oligonucleotide contributes only its best-scoring identification.
        /// </summary>
        /// <seealso cref="BioPolymerGroupExtensions.Score"/>
        public double BioPolymerGroupScore { get; set; }

        /// <summary>
        /// All biopolymer sequences with set modifications identified in this group,
        /// including those shared with other biopolymer groups.
        /// </summary>
        public HashSet<IBioPolymerWithSetMods> AllBioPolymersWithSetMods { get; set; }

        /// <summary>
        /// Biopolymer sequences with set modifications that are unique to this group
        /// (not shared with any other biopolymer group). Used for protein inference.
        /// </summary>
        public HashSet<IBioPolymerWithSetMods> UniqueBioPolymersWithSetMods { get; set; }

        /// <summary>
        /// All peptide-spectrum matches (PSMs) for this group that pass the 1% FDR threshold.
        /// Must be populated before calling <see cref="CalculateSequenceCoverage"/> or <see cref="Score"/>.
        /// Used for scoring, coverage calculation, and quantification.
        /// </summary>
        public HashSet<ISpectralMatch> AllPsmsBelowOnePercentFDR { get; set; }

        /// <summary>
        /// The q-value for this biopolymer group, representing the minimum FDR at which 
        /// this group would be accepted. Lower values indicate higher confidence (0.01 = 1% FDR).
        /// </summary>
        public double QValue { get; set; }

        /// <summary>
        /// The best (lowest) q-value among all biopolymers with set modifications in this group.
        /// </summary>
        public double BestBioPolymerWithSetModsQValue { get; set; }

        /// <summary>
        /// The best (highest) score among all biopolymers with set modifications in this group.
        /// </summary>
        public double BestBioPolymerWithSetModsScore { get; set; }

        /// <summary>
        /// Dictionary mapping sample identifiers to measured intensity values for this group.
        /// Supports both <see cref="SpectraFileInfo"/> (label-free) and <see cref="IsobaricQuantSampleInfo"/> (TMT/iTRAQ) as keys.
        /// </summary>
        public Dictionary<ISampleInfo, double> IntensitiesBySample { get; set; }

        /// <summary>
        /// All biopolymers in this group ordered alphabetically by accession.
        /// Provides a stable, deterministic ordering for output and comparison.
        /// </summary>
        public List<IBioPolymer> ListOfBioPolymersOrderedByAccession { get; private set; }

        #endregion

        #region Additional Properties

        /// <summary>
        /// Cumulative count of target (non-decoy) groups up to and including this one,
        /// when ordered by score descending. Used for FDR calculation via target-decoy approach.
        /// </summary>
        public int CumulativeTarget { get; set; }

        /// <summary>
        /// Cumulative count of decoy groups up to and including this one,
        /// when ordered by score descending. Used for FDR calculation via target-decoy approach.
        /// </summary>
        public int CumulativeDecoy { get; set; }

        /// <summary>
        /// Controls whether modifications are displayed in sequence output.
        /// If true, <see cref="IBioPolymerWithSetMods.FullSequence"/> is used (includes modification annotations).
        /// If false, <see cref="IBioPolymerWithSetMods.BaseSequence"/> is used (unmodified sequence only).
        /// </summary>
        public bool DisplayModsOnPeptides { get; set; }

        /// <summary>
        /// Cached sequence coverage results from <see cref="CalculateSequenceCoverage"/>.
        /// Null until coverage is calculated. Invalidated when <see cref="MergeWith"/> is called.
        /// </summary>
        private SequenceCoverageResult? _coverageResult;

        #endregion

        #region Methods

        /// <summary>
        /// Returns a tab-separated header line for output files, matching the format of <see cref="ToString"/>.
        /// Header includes columns for biopolymer information, quantification, and statistical metrics.
        /// </summary>
        /// <returns>Tab-separated header string suitable for TSV file output.</returns>
        public string GetTabSeparatedHeader()
        {
            var sb = new StringBuilder();
            sb.Append("BioPolymer Accession" + '\t');
            sb.Append("Gene" + '\t');
            sb.Append("Organism" + '\t');
            sb.Append("BioPolymer Full Name" + '\t');
            sb.Append("BioPolymer Unmodified Mass" + '\t');
            sb.Append("Number of BioPolymers in Group" + '\t');
            sb.Append("Unique Sequences" + '\t');
            sb.Append("Shared Sequences" + '\t');
            sb.Append("Number of Sequences" + '\t');
            sb.Append("Number of Unique Sequences" + '\t');
            sb.Append("Sequence Coverage Fraction" + '\t');
            sb.Append("Sequence Coverage" + '\t');
            sb.Append("Sequence Coverage with Mods" + '\t');
            sb.Append("Fragment Sequence Coverage" + '\t');
            sb.Append("Modification Info List" + "\t");

            if (SamplesForQuantification != null)
            {
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
            sb.Append("BioPolymer Decoy/Contaminant/Target" + '\t');
            sb.Append("BioPolymer Cumulative Target" + '\t');
            sb.Append("BioPolymer Cumulative Decoy" + '\t');
            sb.Append("BioPolymer QValue" + '\t');
            sb.Append("Best Sequence Score" + '\t');
            sb.Append("Best Sequence Notch QValue");
            return sb.ToString();
        }

        /// <summary>
        /// Returns a tab-separated string representation of this biopolymer group,
        /// including all identification, quantification, and statistical information.
        /// Format matches the header returned by <see cref="GetTabSeparatedHeader"/>.
        /// </summary>
        /// <returns>Tab-separated string suitable for TSV file output.</returns>
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
            var masses = sequences.Select(sequence =>
                AllBioPolymersWithSetMods.FirstOrDefault(bpws => bpws.BaseSequence == sequence)?.MonoisotopicMass ?? double.NaN);

            sb.Append(TruncateString(string.Join("|", masses)));
            sb.Append("\t");

            sb.Append("" + BioPolymers.Count);
            sb.Append("\t");

            // Compute unique and shared sequences directly
            var (uniqueSeqOutput, sharedSeqOutput) = GetIdentifiedSequencesOutput();
            sb.Append(TruncateString(uniqueSeqOutput));
            sb.Append("\t");
            sb.Append(TruncateString(sharedSeqOutput));
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

            // Use cached coverage results (empty if not yet calculated)
            var coverage = _coverageResult ?? new SequenceCoverageResult();

            sb.Append(TruncateString(string.Join("|",
                coverage.SequenceCoverageFraction.Select(p => string.Format("{0:0.#####}", p)))));
            sb.Append("\t");

            sb.Append(TruncateString(string.Join("|", coverage.SequenceCoverageDisplayList)));
            sb.Append("\t");

            sb.Append(TruncateString(string.Join("|", coverage.SequenceCoverageDisplayListWithMods)));
            sb.Append("\t");

            sb.Append(TruncateString(string.Join("|", coverage.FragmentSequenceCoverageDisplayList)));
            sb.Append("\t");

            sb.Append(TruncateString(string.Join("|", coverage.ModsInfo)));
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

        /// <summary>
        /// Computes pipe-delimited output strings for unique and shared sequences identified in this group.
        /// Uses <see cref="DisplayModsOnPeptides"/> to determine whether to include modification annotations.
        /// </summary>
        /// <returns>Tuple of (uniqueSequences, sharedSequences) pipe-delimited output strings.</returns>
        private (string UniqueSequences, string SharedSequences) GetIdentifiedSequencesOutput()
        {
            var sharedSequences = AllBioPolymersWithSetMods.Except(UniqueBioPolymersWithSetMods);

            string uniqueOutput;
            string sharedOutput;

            if (!DisplayModsOnPeptides)
            {
                uniqueOutput = string.Join("|", UniqueBioPolymersWithSetMods.Select(p => p.BaseSequence).Distinct());
                sharedOutput = string.Join("|", sharedSequences.Select(p => p.BaseSequence).Distinct());
            }
            else
            {
                uniqueOutput = string.Join("|", UniqueBioPolymersWithSetMods.Select(p => p.FullSequence).Distinct());
                sharedOutput = string.Join("|", sharedSequences.Select(p => p.FullSequence).Distinct());
            }

            return (uniqueOutput, sharedOutput);
        }

        /// <summary>
        /// Calculates and updates <see cref="BioPolymerGroupScore"/> based on PSM scores.
        /// Delegates to <see cref="BioPolymerGroupExtensions.Score"/> which computes the score
        /// as the sum of the best (highest) score for each unique base sequence among 
        /// <see cref="AllPsmsBelowOnePercentFDR"/>.
        /// </summary>
        /// <remarks>
        /// This method should be called after <see cref="AllPsmsBelowOnePercentFDR"/> has been populated.
        /// If the collection is empty, <see cref="BioPolymerGroupScore"/> will be set to 0.
        /// </remarks>
        /// <seealso cref="BioPolymerGroupExtensions.Score"/>
        public void Score()
        {
            BioPolymerGroupExtensions.Score(this);
        }

        /// <summary>
        /// Merges another biopolymer group into this one, combining their members, PSMs, and sequences.
        /// Used when groups are determined to represent the same biological entity.
        /// The other group's score is reset to 0 after merging.
        /// </summary>
        /// <param name="otherBioPolymerGroup">The group to merge into this one.</param>
        /// <remarks>
        /// After merging:
        /// <list type="bullet">
        ///   <item><description><see cref="BioPolymers"/> contains the union of both groups' biopolymers</description></item>
        ///   <item><description><see cref="AllBioPolymersWithSetMods"/> contains the union of both groups' sequences</description></item>
        ///   <item><description><see cref="UniqueBioPolymersWithSetMods"/> contains the union of both groups' unique sequences</description></item>
        ///   <item><description><see cref="AllPsmsBelowOnePercentFDR"/> contains the union of both groups' PSMs</description></item>
        ///   <item><description><see cref="ListOfBioPolymersOrderedByAccession"/> and <see cref="BioPolymerGroupName"/> are recalculated</description></item>
        ///   <item><description>Cached coverage results are invalidated</description></item>
        /// </list>
        /// </remarks>
        public void MergeWith(IBioPolymerGroup otherBioPolymerGroup)
        {
            this.BioPolymers.UnionWith(otherBioPolymerGroup.BioPolymers);
            this.AllBioPolymersWithSetMods.UnionWith(otherBioPolymerGroup.AllBioPolymersWithSetMods);
            this.UniqueBioPolymersWithSetMods.UnionWith(otherBioPolymerGroup.UniqueBioPolymersWithSetMods);
            this.AllPsmsBelowOnePercentFDR.UnionWith(otherBioPolymerGroup.AllPsmsBelowOnePercentFDR);
            otherBioPolymerGroup.BioPolymerGroupScore = 0;

            ListOfBioPolymersOrderedByAccession = BioPolymers.OrderBy(p => p.Accession).ToList();
            BioPolymerGroupName = string.Join("|", ListOfBioPolymersOrderedByAccession.Select(p => p.Accession));

            // Invalidate cached coverage since PSMs changed
            _coverageResult = null;
        }

        /// <summary>
        /// Creates a new biopolymer group containing only data from a specific spectra file.
        /// The subset group will have the same biopolymers but filtered PSMs, sequences, and intensities.
        /// Used for per-file analysis and output.
        /// </summary>
        /// <param name="fullFilePath">The full path to the spectra file to filter by.</param>
        /// <param name="silacLabels">Optional SILAC labels to apply during subset creation.</param>
        /// <returns>A new <see cref="BioPolymerGroup"/> containing only data from the specified file.</returns>
        public IBioPolymerGroup ConstructSubsetBioPolymerGroup(string fullFilePath, List<SilacLabel>? silacLabels = null)
        {
            var allPsmsForThisFile =
                new HashSet<ISpectralMatch>(
                    AllPsmsBelowOnePercentFDR.Where(p => p.FullFilePath.Equals(fullFilePath)));
            var allSequencesForThisFile =
                new HashSet<IBioPolymerWithSetMods>(
                    allPsmsForThisFile.SelectMany(p => p.GetIdentifiedBioPolymersWithSetMods()));
            var allUniqueSequencesForThisFile =
                new HashSet<IBioPolymerWithSetMods>(UniqueBioPolymersWithSetMods.Intersect(allSequencesForThisFile));

            BioPolymerGroup subsetGroup = new BioPolymerGroup(BioPolymers, allSequencesForThisFile, allUniqueSequencesForThisFile)
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

                subsetGroup.SamplesForQuantification = matchingSamples;

                if (IntensitiesBySample != null)
                {
                    subsetGroup.IntensitiesBySample = IntensitiesBySample
                        .Where(kvp => matchingSamples.Contains(kvp.Key))
                        .ToDictionary(kvp => kvp.Key, kvp => kvp.Value);
                }
            }

            return subsetGroup;
        }

        /// <summary>
        /// Calculates sequence coverage for all biopolymers in this group at two levels:
        /// <list type="number">
        ///   <item><description><b>Peptide-level coverage:</b> All residues within identified peptide boundaries 
        ///   are considered covered. Results cached in <see cref="_coverageResult"/>.</description></item>
        ///   <item><description><b>Fragment-level coverage:</b> Only residues with supporting fragment ion evidence 
        ///   are considered covered. Requires PSMs to implement <see cref="IHasSequenceCoverageFromFragments"/>.</description></item>
        /// </list>
        /// 
        /// Display strings use uppercase letters for covered residues and lowercase for uncovered residues.
        /// Also calculates modification occupancy statistics and populates <see cref="SequenceCoverageResult.ModsInfo"/>.
        /// </summary>
        /// <remarks>
        /// This method should be called after <see cref="AllPsmsBelowOnePercentFDR"/> has been populated.
        /// If PSMs do not implement <see cref="IHasSequenceCoverageFromFragments"/>, fragment-level coverage
        /// will show all residues as uncovered (all lowercase).
        /// Results are cached and used by <see cref="ToString"/> for output formatting.
        /// </remarks>
        public void CalculateSequenceCoverage()
        {
            var result = new SequenceCoverageResult();

            // Maps biopolymers to their identified sequences with unambiguous base sequences
            var bioPolymersWithUnambiguousSequences = new Dictionary<IBioPolymer, List<IBioPolymerWithSetMods>>();
            // Maps biopolymers to sequences with successfully localized modifications
            var bioPolymersWithLocalizedMods = new Dictionary<IBioPolymer, List<IBioPolymerWithSetMods>>();

            foreach (var bioPolymer in BioPolymers)
            {
                bioPolymersWithUnambiguousSequences.Add(bioPolymer, new List<IBioPolymerWithSetMods>());
                bioPolymersWithLocalizedMods.Add(bioPolymer, new List<IBioPolymerWithSetMods>());
            }

            // Check once if PSMs support fragment coverage calculation
            bool supportsFragmentCoverage = AllPsmsBelowOnePercentFDR.FirstOrDefault() is IHasSequenceCoverageFromFragments;

            // If fragment coverage is supported, calculate it for all PSMs upfront
            if (supportsFragmentCoverage)
            {
                foreach (var psm in AllPsmsBelowOnePercentFDR.Where(p => p.BaseSequence != null))
                {
                    ((IHasSequenceCoverageFromFragments)psm).GetSequenceCoverage();
                }
            }

            // Collect sequences from PSMs with unambiguous identifications
            foreach (var psm in AllPsmsBelowOnePercentFDR.Where(p => p.BaseSequence != null))
            {
                foreach (var sequence in psm.GetIdentifiedBioPolymersWithSetMods()
                    .DistinctBy(p => p.FullSequence)
                    .Where(s => BioPolymers.Contains(s.Parent)))
                {
                    bioPolymersWithUnambiguousSequences[sequence.Parent].Add(sequence);

                    // null FullSequence means mods were not localized; don't include in mods display
                    if (sequence.FullSequence != null)
                    {
                        bioPolymersWithLocalizedMods[sequence.Parent].Add(sequence);
                    }
                }
            }

            // Calculate fragment-level sequence coverage (amino acid level based on fragment ions)
            foreach (var bioPolymer in ListOfBioPolymersOrderedByAccession)
            {
                var coveredResiduesOneBased = new HashSet<int>();

                // Only process fragment coverage if PSMs support it
                if (supportsFragmentCoverage)
                {
                    foreach (var psm in AllPsmsBelowOnePercentFDR.Where(p => p.BaseSequence != null))
                    {
                        var coverageProvider = (IHasSequenceCoverageFromFragments)psm;

                        if (coverageProvider.FragmentCoveragePositionInPeptide == null)
                            continue;

                        // Get sequences from this PSM that belong to this biopolymer
                        var sequencesForThisBioPolymer = psm.GetIdentifiedBioPolymersWithSetMods()
                            .Where(p => p.Parent.Accession == bioPolymer.Accession);

                        foreach (var sequence in sequencesForThisBioPolymer)
                        {
                            // Convert peptide positions to protein positions
                            foreach (var position in coverageProvider.FragmentCoveragePositionInPeptide)
                            {
                                // Convert a one-based peptide position to a one-based protein position
                                // by adding the peptide's starting residue in the protein and subtracting 1.
                                // This accounts for the peptide's offset within the protein sequence.
                                int proteinPosition = position + sequence.OneBasedStartResidue - 1;
                                coveredResiduesOneBased.Add(proteinPosition);
                            }
                        }
                    }
                }

                // Build display string: uppercase = covered, lowercase = not covered
                char[] fragmentCoverageArray = bioPolymer.BaseSequence.ToLower().ToCharArray();
                foreach (var residue in coveredResiduesOneBased.Where(r => r >= 1 && r <= fragmentCoverageArray.Length))
                {
                    fragmentCoverageArray[residue - 1] = char.ToUpper(fragmentCoverageArray[residue - 1]);
                }

                result.FragmentSequenceCoverageDisplayList.Add(new string(fragmentCoverageArray));
            }

            // Calculate peptide-level sequence coverage (all residues in identified peptides are covered)
            foreach (var bioPolymer in ListOfBioPolymersOrderedByAccession)
            {
                var coveredResiduesOneBased = new HashSet<int>();

                // Mark all residues within each identified peptide as covered
                foreach (var sequence in bioPolymersWithUnambiguousSequences[bioPolymer])
                {
                    for (int i = sequence.OneBasedStartResidue; i <= sequence.OneBasedEndResidue; i++)
                    {
                        coveredResiduesOneBased.Add(i);
                    }
                }

                // Calculate coverage fraction
                double coverageFraction = (double)coveredResiduesOneBased.Count / bioPolymer.Length;
                result.SequenceCoverageFraction.Add(coverageFraction);

                // Build display string: uppercase = covered, lowercase = not covered
                char[] coverageArray = bioPolymer.BaseSequence.ToLower().ToCharArray();
                foreach (var residueLocation in coveredResiduesOneBased.Where(r => r >= 1 && r <= coverageArray.Length))
                {
                    coverageArray[residueLocation - 1] = char.ToUpper(coverageArray[residueLocation - 1]);
                }

                string sequenceCoverageDisplay = new string(coverageArray);
                result.SequenceCoverageDisplayList.Add(sequenceCoverageDisplay);

                // Build coverage display with modifications
                var modsOnThisBioPolymer = new HashSet<KeyValuePair<int, Modification>>();

                foreach (var sequence in bioPolymersWithLocalizedMods[bioPolymer])
                {
                    foreach (var mod in sequence.AllModsOneIsNterminus)
                    {
                        // Skip peptide terminal mods and common variable/fixed mods
                        if (mod.Value.ModificationType.Contains("PeptideTermMod") ||
                            mod.Value.ModificationType.Contains("Common Variable") ||
                            mod.Value.ModificationType.Contains("Common Fixed"))
                        {
                            continue;
                        }

                        // Convert from AllModsOneIsNterminus indexing (where 1 is N-terminus) to protein position.
                        // Subtracting 2 aligns the peptide-local modification position to the protein sequence.
                        int proteinPosition = sequence.OneBasedStartResidue + mod.Key - 2;
                        modsOnThisBioPolymer.Add(new KeyValuePair<int, Modification>(proteinPosition, mod.Value));
                    }
                }

                // Insert modification annotations into sequence coverage display
                var sequenceCoverageWithModsBuilder = new StringBuilder(sequenceCoverageDisplay);
                var orderedMods = modsOnThisBioPolymer.OrderBy(p => p.Key).ToList();

                // Track offset as we insert modifications (each insertion shifts subsequent positions)
                int insertionOffset = 0;

                foreach (var mod in orderedMods)
                {
                    if (mod.Value.LocationRestriction.Equals("N-terminal."))
                    {
                        string prefix = $"[{mod.Value.IdWithMotif}]-";
                        sequenceCoverageWithModsBuilder.Insert(0, prefix);
                        insertionOffset += prefix.Length;
                    }
                    else if (mod.Value.LocationRestriction.Equals("Anywhere."))
                    {
                        int baseInsertIndex = sequenceCoverageDisplay.Length - (bioPolymer.Length - mod.Key);
                        if (baseInsertIndex >= 0 && baseInsertIndex <= sequenceCoverageDisplay.Length)
                        {
                            string modAnnotation = $"[{mod.Value.IdWithMotif}]";
                            sequenceCoverageWithModsBuilder.Insert(baseInsertIndex + insertionOffset, modAnnotation);
                            insertionOffset += modAnnotation.Length;
                        }
                    }
                    else if (mod.Value.LocationRestriction.Equals("C-terminal."))
                    {
                        sequenceCoverageWithModsBuilder.Append($"-[{mod.Value.IdWithMotif}]");
                    }
                }

                result.SequenceCoverageDisplayListWithMods.Add(sequenceCoverageWithModsBuilder.ToString());

                // Calculate modification occupancy statistics
                if (modsOnThisBioPolymer.Any())
                {
                    CalculateModificationOccupancy(bioPolymer, bioPolymersWithLocalizedMods[bioPolymer], result);
                }
            }

            _coverageResult = result;
        }

        /// <summary>
        /// Calculates modification occupancy statistics for a biopolymer and adds them to <see cref="SequenceCoverageResult.ModsInfo"/>.
        /// Occupancy is calculated as the fraction of peptides covering each modification site that contain the modification.
        /// </summary>
        /// <param name="bioPolymer">The biopolymer to calculate modification occupancy for.</param>
        /// <param name="localizedSequences">Sequences with localized modifications mapping to this biopolymer.</param>
        /// <param name="result">The result object to populate with modification info.</param>
        private void CalculateModificationOccupancy(IBioPolymer bioPolymer, List<IBioPolymerWithSetMods> localizedSequences, SequenceCoverageResult result)
        {
            var modCounts = new List<int>();      // Count of modified peptides at each position
            var totalCounts = new List<int>();    // Count of all peptides covering each position
            var modPositions = new List<(int index, string modName)>();

            foreach (var sequence in localizedSequences)
            {
                foreach (var mod in sequence.AllModsOneIsNterminus)
                {
                    // Skip common variable/fixed mods and peptide terminal mods
                    if (mod.Value.ModificationType.Contains("Common Variable") ||
                        mod.Value.ModificationType.Contains("Common Fixed") ||
                        mod.Value.LocationRestriction.Equals("NPep") ||
                        mod.Value.LocationRestriction.Equals("PepC"))
                    {
                        continue;
                    }

                    int indexInProtein;
                    if (mod.Value.LocationRestriction.Equals("N-terminal."))
                    {
                        indexInProtein = 1;
                    }
                    else if (mod.Value.LocationRestriction.Equals("Anywhere."))
                    {
                        indexInProtein = sequence.OneBasedStartResidue + mod.Key - 2;
                    }
                    else if (mod.Value.LocationRestriction.Equals("C-terminal."))
                    {
                        indexInProtein = bioPolymer.Length;
                    }
                    else
                    {
                        // Skip unrecognized location restrictions
                        continue;
                    }

                    var modKey = (indexInProtein, mod.Value.IdWithMotif);

                    if (modPositions.Contains(modKey))
                    {
                        modCounts[modPositions.IndexOf(modKey)]++;
                    }
                    else
                    {
                        modPositions.Add(modKey);

                        // Count total peptides covering this position
                        int peptidesAtPosition = 0;
                        foreach (var seq in localizedSequences)
                        {
                            int rangeStart = seq.OneBasedStartResidue - (indexInProtein == 1 ? 1 : 0);
                            if (indexInProtein >= rangeStart && indexInProtein <= seq.OneBasedEndResidue)
                            {
                                peptidesAtPosition++;
                            }
                        }

                        totalCounts.Add(peptidesAtPosition);
                        modCounts.Add(1);
                    }
                }
            }

            // Build modification info string
            var modStrings = new List<(int position, string info)>();
            for (int i = 0; i < modCounts.Count; i++)
            {
                string position = modPositions[i].index.ToString();
                string modName = modPositions[i].modName;
                string occupancy = ((double)modCounts[i] / totalCounts[i]).ToString("F2");
                string fractionalOccupancy = $"{modCounts[i]}/{totalCounts[i]}";
                string modString = $"#aa{position}[{modName},info:occupancy={occupancy}({fractionalOccupancy})]";
                modStrings.Add((modPositions[i].index, modString));
            }

            string modInfoString = string.Join(";", modStrings.OrderBy(x => x.position).Select(x => x.info));

            if (!string.IsNullOrEmpty(modInfoString))
            {
                result.ModsInfo.Add(modInfoString);
            }
        }

        #endregion

        #region Equality

        /// <summary>
        /// Determines whether this biopolymer group equals another based on group name.
        /// Two groups are considered equal if they have the same <see cref="BioPolymerGroupName"/>.
        /// </summary>
        /// <param name="other">The other biopolymer group to compare.</param>
        /// <returns>True if the groups have the same name; otherwise, false.</returns>
        public bool Equals(IBioPolymerGroup? other)
        {
            if (other is null) return false;
            if (ReferenceEquals(this, other)) return true;
            return BioPolymerGroupName == other.BioPolymerGroupName;
        }

        /// <summary>
        /// Determines whether this biopolymer group equals another object.
        /// Supports comparison with both <see cref="IBioPolymerGroup"/> and <see cref="BioPolymerGroup"/> types.
        /// </summary>
        /// <param name="obj">The object to compare.</param>
        /// <returns>True if the objects are equal; otherwise, false.</returns>
        public override bool Equals(object? obj)
        {
            if (obj is BioPolymerGroup bg) return Equals(bg);
            if (obj is IBioPolymerGroup ibg) return Equals(ibg);
            return false;
        }

        /// <summary>
        /// Returns a hash code for this biopolymer group based on <see cref="BioPolymerGroupName"/>.
        /// </summary>
        /// <returns>Hash code integer.</returns>
        public override int GetHashCode()
        {
            return BioPolymerGroupName?.GetHashCode() ?? 0;
        }

        #endregion

        #region Private Helpers

        /// <summary>
        /// Truncates a string to <see cref="MaxStringLength"/> if it exceeds that length.
        /// Used to ensure output compatibility with Excel's cell character limits.
        /// </summary>
        /// <param name="input">The string to truncate.</param>
        /// <returns>Truncated string, original string if within limits, or empty string if input is null/empty.</returns>
        private static string TruncateString(string? input)
        {
            if (string.IsNullOrEmpty(input))
                return string.Empty;

            if (MaxStringLength <= 0 || input.Length <= MaxStringLength)
                return input;

            return input.Substring(0, MaxStringLength);
        }

        /// <summary>
        /// Holds cached sequence coverage calculation results from <see cref="CalculateSequenceCoverage"/>.
        /// Encapsulates the various coverage display lists to avoid storing them as separate class properties.
        /// </summary>
        public sealed class SequenceCoverageResult
        {
            /// <summary>
            /// Sequence coverage fraction for each biopolymer in the group, ordered by accession.
            /// Each value (0.0 to 1.0) represents the fraction of residues covered by identified peptides.
            /// </summary>
            public List<double> SequenceCoverageFraction { get; } = new();

            /// <summary>
            /// Visual representation of sequence coverage for each biopolymer in the group, ordered by accession.
            /// Uppercase letters indicate covered residues; lowercase indicates uncovered residues.
            /// </summary>
            public List<string> SequenceCoverageDisplayList { get; } = new();

            /// <summary>
            /// Visual representation of sequence coverage including modification annotations, ordered by accession.
            /// Modifications are shown as [ModName] inserted at the appropriate position.
            /// </summary>
            public List<string> SequenceCoverageDisplayListWithMods { get; } = new();

            /// <summary>
            /// Visual representation of fragment-level sequence coverage for each biopolymer, ordered by accession.
            /// Uppercase letters indicate residues covered by matched fragment ions; lowercase indicates uncovered.
            /// Will show all lowercase if PSMs do not implement <see cref="IHasSequenceCoverageFromFragments"/>.
            /// </summary>
            public List<string> FragmentSequenceCoverageDisplayList { get; } = new();

            /// <summary>
            /// Modification occupancy information for this group. Each string describes a modification 
            /// at a specific position with its occupancy fraction (e.g., how often the site is modified).
            /// Format: #aa{position}[{modName},info:occupancy={fraction}({count}/{total})]
            /// </summary>
            public List<string> ModsInfo { get; } = new();
        }

        #endregion
    }
}