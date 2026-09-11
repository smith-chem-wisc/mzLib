using Easy.Common.Extensions;
using MassSpectrometry;
using MzLibUtil;
using Omics.Modifications;
using Omics.SpectralMatch;
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
    /// </list>
    ///
    /// TSV output is defined by <see cref="BioPolymerGroupTsvSchema"/> and rendered by
    /// <see cref="TsvWriter"/>. The compatibility members <see cref="GetTabSeparatedHeader"/> and
    /// <see cref="ToString"/> delegate to a schema built from this group alone. A writer rendering a
    /// whole file should instead build one schema from every group, because the quantification
    /// columns depend on the complete dataset and every row must use the same schema.
    ///
    /// Fragment-level coverage is only calculated when PSMs implement <see cref="IHasSequenceCoverageFromFragments"/>.
    /// </summary>
    public class BioPolymerGroup : IBioPolymerGroup
    {
        /// <summary>
        /// Creates a new biopolymer group from the specified biopolymers and identified sequences.
        /// </summary>
        /// <param name="bioPolymers">Set of biopolymers (e.g., proteins, RNA) that belong to this group.
        /// These are typically indistinguishable based on the identified sequences.</param>
        /// <param name="bioPolymersWithSetMods">All identified sequences with modifications for this group,
        /// including sequences shared with other groups.</param>
        /// <param name="uniqueBioPolymersWithSetMods">Sequences with modifications that are unique to this group
        /// and not shared with any other biopolymer group.</param>
        /// <param name="groupType">Identifies the type of biopolymer in this group, which determines the modification
        /// occupancy calculation strategy used by <see cref="PopulateSampleGroupResults"/>.
        /// <see cref="BioPolymerGroupType.Parent"/> uses parent(typically protein)-level coordinates;
        /// <see cref="BioPolymerGroupType.DigestionProduct"/> uses
        /// digestion-product-local coordinates (typically peptide positions).</param>
        public BioPolymerGroup(HashSet<IBioPolymer> bioPolymers, HashSet<IBioPolymerWithSetMods> bioPolymersWithSetMods,
            HashSet<IBioPolymerWithSetMods> uniqueBioPolymersWithSetMods, BioPolymerGroupType groupType = BioPolymerGroupType.Parent)
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
            IsEntrapment = false;
            GroupType = groupType;

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

                if (bioPolymer.IsEntrapment)
                {
                    IsEntrapment = true;
                }

                // If all three are true, we can break early
                if (IsDecoy && IsContaminant && IsEntrapment)
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
        /// True if this group contains any biopolymers marked as entrapment proteins.
        /// </summary>
        public bool IsEntrapment { get; }

        /// <summary>
        /// List of samples that contribute quantification data for this group.
        /// Supports both <see cref="SpectraFileInfo"/> (label-free) and <see cref="IsobaricQuantSampleInfo"/> (TMT/iTRAQ).
        /// Setting this property invalidates <see cref="SampleGroupResults"/>, which will be 
        /// re-populated on the next call to <see cref="PopulateSampleGroupResults"/>.
        /// </summary>
        private List<ISampleInfo>? _samplesForQuantification;
        public List<ISampleInfo>? SamplesForQuantification
        {
            get => _samplesForQuantification;
            set
            {
                _samplesForQuantification = value;
                SampleGroupResults = null;
            }
        }

        /// <summary>
        /// Dictionary mapping sample identifiers to measured intensity values for this group.
        /// Supports both <see cref="SpectraFileInfo"/> (label-free) and <see cref="IsobaricQuantSampleInfo"/> (TMT/iTRAQ) as keys.
        /// Setting this property invalidates <see cref="SampleGroupResults"/>, which will be
        /// re-populated on the next call to <see cref="PopulateSampleGroupResults"/>.
        /// </summary>
        private Dictionary<ISampleInfo, double>? _intensitiesBySample;
        public Dictionary<ISampleInfo, double>? IntensitiesBySample
        {
            get => _intensitiesBySample;
            set
            {
                _intensitiesBySample = value;
                SampleGroupResults = null;
            }
        }

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
        /// Aggregated confidence score for the group, used internally for protein grouping optimization.
        /// Computed by <see cref="Score"/> as the sum of the best (highest) score for each unique 
        /// base sequence among the PSMs in <see cref="AllPsmsBelowOnePercentFDR"/>.
        /// Higher values indicate higher confidence. NOT used for protein FDR calculations.
        /// </summary>
        /// <seealso cref="Score"/>
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
        /// Setting this property invalidates both <see cref="SampleGroupResults"/> and the cached
        /// sequence coverage result.
        /// </summary>
        private HashSet<ISpectralMatch> _allPsmsBelowOnePercentFDR = null!;
        public HashSet<ISpectralMatch> AllPsmsBelowOnePercentFDR
        {
            get => _allPsmsBelowOnePercentFDR;
            set
            {
                _allPsmsBelowOnePercentFDR = value;
                SampleGroupResults = null;
                _coverageResult = null;
            }
        }

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
        /// All biopolymers in this group ordered alphabetically by accession.
        /// Provides a stable, deterministic ordering for output and comparison.
        /// </summary>
        public List<IBioPolymer> ListOfBioPolymersOrderedByAccession { get; private set; }

        /// <summary>
        /// Per-sample-group quantification and modification occupancy results.
        /// Each entry represents one (Condition × BiologicalReplicate) group for label-free data,
        /// or one (File × Channel) for isobaric data.
        /// Built by <see cref="PopulateSampleGroupResults"/> from <see cref="SamplesForQuantification"/>,
        /// <see cref="IntensitiesBySample"/>, and <see cref="AllPsmsBelowOnePercentFDR"/>.
        /// Consumed by <see cref="BioPolymerGroupTsvSchema"/> to build the per-sample output columns.
        /// </summary>
        public List<SampleGroupResult>? SampleGroupResults { get; set; }

        #endregion

        #region Rendering

        /// <summary>
        /// Tab-separated header line describing this group's own columns, matching the row returned
        /// by <see cref="ToString"/>.
        ///
        /// Kept on the type, delegating to <see cref="BioPolymerGroupTsvSchema"/>, so there is one
        /// implementation of how a column is produced while the member stays where its callers
        /// expect it. MetaMorpheus's ProteinGroup overrides this and reads
        /// <see cref="HasAssignedSampleIntensities"/>, so removing either in the same release that
        /// introduces the schema would make that release one MetaMorpheus cannot adopt.
        ///
        /// A writer rendering a whole file should build the schema once from every group
        /// (<see cref="BioPolymerGroupTsvSchema.For"/>) rather than taking the header from one group
        /// and rows from the rest -- deriving header and rows from separate groups is what made the
        /// table ragged before, and the dataset-wide schema is what fixes it.
        /// </summary>
        public virtual string GetTabSeparatedHeader() =>
            TsvWriter.HeaderLine(BioPolymerGroupTsvSchema.For(new[] { this }));

        /// <summary>
        /// This group rendered as one tab-separated row, aligned with
        /// <see cref="GetTabSeparatedHeader"/>.
        /// </summary>
        public override string ToString() =>
            TsvWriter.RowLine(BioPolymerGroupTsvSchema.For(new[] { this }), this);

        /// <summary>
        /// Whether an engine assigned quantification to this group, and therefore whether its
        /// intensity columns exist at all.
        ///
        /// Deliberately NOT <see cref="SampleGroupResult.HasIntensityData"/>. That answers a
        /// per-sample-group question -- "did this sample group receive a value" -- and using it to
        /// choose columns made a group whose samples were all unobserved describe a narrower table
        /// than its neighbour. Both fields are required: an intensity dictionary without samples has
        /// nothing to label, and samples without a dictionary have nothing to report.
        ///
        /// This is the same predicate <see cref="BioPolymerGroupTsvSchema.For"/> applies across a
        /// dataset, evaluated over this group alone -- so the single-group schema above reproduces
        /// exactly what this member used to gate, and no rendered output changes.
        /// </summary>
        protected bool HasAssignedSampleIntensities =>
            SamplesForQuantification is { Count: > 0 } && IntensitiesBySample is not null;

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
        /// Identifies the type of biopolymer in this group, which determines the modification
        /// occupancy calculation strategy used by <see cref="PopulateSampleGroupResults"/>.
        /// <see cref="BioPolymerGroupType.Parent"/> uses protein-level coordinates;
        /// <see cref="BioPolymerGroupType.DigestionProduct"/> use digestion-product-local coordinates.
        /// </summary>
        public BioPolymerGroupType GroupType { get; }

        /// <summary>
        /// Cached sequence coverage results from <see cref="CalculateSequenceCoverage"/>.
        /// Null until coverage is calculated. Invalidated when <see cref="MergeWith"/> is called.
        /// </summary>
        private SequenceCoverageResult? _coverageResult;
        public SequenceCoverageResult CoverageResult
        {
            get
            {
                if (_coverageResult is null)
                    CalculateSequenceCoverage();
                return _coverageResult!;
            }
        }

        /// <summary>
        /// True once <see cref="CalculateSequenceCoverage"/> has run and its result is still valid.
        /// Lets a caller read <see cref="CoverageResult"/> only when it is already available, since
        /// reading it otherwise triggers the calculation.
        /// </summary>
        public bool IsSequenceCoverageCalculated => _coverageResult is not null;

        #endregion

        #region Methods

        /// <summary>
        /// Builds <see cref="SampleGroupResults"/> from the existing <see cref="SamplesForQuantification"/>,
        /// <see cref="IntensitiesBySample"/>, and <see cref="AllPsmsBelowOnePercentFDR"/>.
        /// Groups samples by (Condition, BiologicalReplicate) for label-free data, by
        /// (File, Channel) for isobaric data, or by PSM file path when no experimental design is available.
        /// For each group, computes spectral count, per-file intensities (stored on the result),
        /// and modification occupancy at both protein and peptide levels.
        /// </summary>
        /// <remarks>
        /// Must be called after <see cref="AllPsmsBelowOnePercentFDR"/> has been populated.
        /// Invoked by <see cref="BioPolymerGroupTsvSchema"/> whenever <see cref="SampleGroupResults"/>
        /// is null — which occurs after construction or after setting
        /// <see cref="SamplesForQuantification"/>, <see cref="IntensitiesBySample"/>,
        /// or <see cref="AllPsmsBelowOnePercentFDR"/>.
        /// </remarks>
        public void PopulateSampleGroupResults()
        {
            SampleGroupResults = SampleGroupBuilder.Build(
                SamplesForQuantification,
                IntensitiesBySample,
                AllPsmsBelowOnePercentFDR,
                PopulateOccupancy);
        }

        /// <summary>
        /// Populates protein-level and peptide-level modification occupancy on a <see cref="SampleGroupResult"/>
        /// using the specified PSMs. PSM grouping, form filtering, TotalCount derivation, and intensity
        /// lookup are all handled internally by <see cref="ModificationOccupancyCalculator"/>.
        /// </summary>
        private void PopulateOccupancy(SampleGroupResult result, List<ISpectralMatch> psms)
        {
            if (GroupType == BioPolymerGroupType.Parent)
            {
                foreach (var bioPolymer in ListOfBioPolymersOrderedByAccession)
                {
                    var occupancy = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
                        bioPolymer, psms);

                    if (occupancy.Count > 0)
                        result.ParentOccupancy[bioPolymer.Accession] = occupancy;
                }
            }
            else
            {
                var psmsGroupedByBaseSequence = psms.GroupBy(p => p.BaseSequence);
                foreach (var baseSeqGroup in psmsGroupedByBaseSequence)
                { 
                    var occupancy = ModificationOccupancyCalculator.CalculateDigestionProductLevelOccupancy(baseSeqGroup.ToList());

                    if (occupancy.Count > 0)
                    {
                        result.DigestionProductOccupancy[baseSeqGroup.Key] = occupancy;
                    }
                }
            }
        }

        /// <summary>
        /// Calculates and updates <see cref="BioPolymerGroupScore"/> based on PSM scores.
        /// 
        /// The score is computed as the sum of the best (highest) score for each unique base sequence
        /// among <see cref="AllPsmsBelowOnePercentFDR"/>. This ensures each unique peptide/oligonucleotide
        /// contributes only its best-scoring identification to the group score.
        /// </summary>
        /// <remarks>
        /// This method is used internally for protein grouping optimization and is NOT used for 
        /// protein FDR calculations.
        /// 
        /// This method should be called after <see cref="AllPsmsBelowOnePercentFDR"/> has been populated.
        /// If the collection is empty, <see cref="BioPolymerGroupScore"/> will be set to 0.
        /// </remarks>
        public void Score()
        {
            BioPolymerGroupScore = AllPsmsBelowOnePercentFDR
                .GroupBy(p => p.BaseSequence)
                .Select(p => p.Select(x => x.Score).Max())
                .Sum();
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
            SampleGroupResults = null;
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

            // ConstructSubsetBioPolymerGroup passes it through the constructor instead of object initializer
            BioPolymerGroup subsetGroup = new BioPolymerGroup(
                BioPolymers,
                allSequencesForThisFile,
                allUniqueSequencesForThisFile,
                GroupType)
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
        /// Display strings use uppercase letters for covered residues and lowercase for uncovered residues.
        /// </summary>
        /// <remarks>
        /// Must be called after <see cref="AllPsmsBelowOnePercentFDR"/> has been populated.
        /// If PSMs do not implement <see cref="IHasSequenceCoverageFromFragments"/>, fragment-level
        /// coverage will show all residues as uncovered (all lowercase).
        /// Results are cached and invalidated by <see cref="MergeWith"/> or reassignment of
        /// <see cref="AllPsmsBelowOnePercentFDR"/>.
        /// </remarks>
        public virtual void CalculateSequenceCoverage()
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
            }

            _coverageResult = result;
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
        }

        #endregion
    }
}
