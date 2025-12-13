using Easy.Common.Extensions;
using MassSpectrometry;
using Omics.Modifications;
using System.Text;

namespace Omics.BioPolymerGroup
{
	/// <summary>
	/// Represents a group of related biopolymers (e.g., proteins, RNA sequences) that share identified sequences/fragments.
	/// Provides quantification, scoring, and output formatting for the group.
	/// This is a generic implementation of IBioPolymerGroup suitable for any biopolymer type.
	/// </summary>
	public class BioPolymerGroup : IBioPolymerGroup
	{
		/// <summary>
		/// Maximum length for string fields in output. Strings longer than this will be truncated.
		/// Set to 0 or negative to disable truncation.
		/// </summary>
		public static int MaxStringLength { get; set; } = 32000;

		/// <summary>
		/// Creates a new biopolymer group from the specified biopolymers and identified sequences.
		/// </summary>
		/// <param name="bioPolymers">Set of biopolymers (e.g., proteins, RNA) that belong to this group.</param>
		/// <param name="bioPolymersWithSetMods">All identified sequences with modifications for this group.</param>
		/// <param name="uniqueBioPolymersWithSetMods">Sequences with modifications unique to this group (not shared with other groups).</param>
		public BioPolymerGroup(HashSet<IBioPolymer> bioPolymers, HashSet<IBioPolymerWithSetMods> bioPolymersWithSetMods,
			HashSet<IBioPolymerWithSetMods> uniqueBioPolymersWithSetMods)
		{
			BioPolymers = bioPolymers;
			ListOfBioPolymersOrderedByAccession = BioPolymers.OrderBy(p => p.Accession).ToList();
			BioPolymerGroupName = string.Join("|", ListOfBioPolymersOrderedByAccession.Select(p => p.Accession));
			AllBioPolymersWithSetMods = bioPolymersWithSetMods;
			UniqueBioPolymersWithSetMods = uniqueBioPolymersWithSetMods;
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
		/// True if this group contains only decoy biopolymers, used for FDR estimation.
		/// </summary>
		public bool IsDecoy { get; }

		/// <summary>
		/// True if this group contains biopolymers marked as contaminants.
		/// </summary>
		public bool IsContaminant { get; }

		/// <summary>
		/// List of samples that contribute quantification data for this group.
		/// Supports both SpectraFileInfo (label-free) and IsobaricQuantSampleInfo (TMT/iTRAQ).
		/// </summary>
		public List<ISampleInfo> SamplesForQuantification { get; set; }

		/// <summary>
		/// Set of all biopolymers (e.g., proteins, RNA sequences) that belong to this group.
		/// </summary>
		public HashSet<IBioPolymer> BioPolymers { get; set; }

		/// <summary>
		/// Display name for the biopolymer group, derived from the accessions of member biopolymers.
		/// Used as the primary identity key for equality comparisons.
		/// </summary>
		public string BioPolymerGroupName { get; private set; }

		/// <summary>
		/// Aggregated score for the group, computed from member PSMs.
		/// Higher scores typically indicate higher confidence in the identification.
		/// </summary>
		public double BioPolymerGroupScore { get; set; }

		/// <summary>
		/// All biopolymer sequences with set modifications identified in this group,
		/// including those shared with other groups.
		/// </summary>
		public HashSet<IBioPolymerWithSetMods> AllBioPolymersWithSetMods { get; set; }

		/// <summary>
		/// Biopolymer sequences with set modifications that are unique to this group
		/// (not shared with any other biopolymer group).
		/// </summary>
		public HashSet<IBioPolymerWithSetMods> UniqueBioPolymersWithSetMods { get; set; }

		/// <summary>
		/// All peptide-spectrum matches (PSMs) for this group that pass the 1% FDR threshold.
		/// </summary>
		public HashSet<ISpectralMatch> AllPsmsBelowOnePercentFDR { get; set; }

		/// <summary>
		/// The q-value (FDR-adjusted p-value) for this biopolymer group.
		/// Lower values indicate higher confidence in the identification.
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
		/// Summary information about modifications present on members of this group.
		/// Each string typically describes a modification type and its frequency or location.
		/// </summary>
		public List<string> ModsInfo { get; private set; }

		/// <summary>
		/// Dictionary mapping sample identifiers to measured intensity values for this group.
		/// Supports both SpectraFileInfo (label-free) and IsobaricQuantSampleInfo (TMT/iTRAQ) as keys.
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
		/// Sequence coverage fraction for each biopolymer in the group, ordered by accession.
		/// Each value represents the fraction of the biopolymer sequence covered by identified fragments.
		/// </summary>
		public List<double> SequenceCoverageFraction { get; private set; }

		/// <summary>
		/// Visual representation of sequence coverage for each biopolymer in the group.
		/// Typically shows identified regions with special characters or formatting.
		/// </summary>
		public List<string> SequenceCoverageDisplayList { get; private set; }

		/// <summary>
		/// Visual representation of sequence coverage including modification information.
		/// Shows both coverage and locations of identified modifications.
		/// </summary>
		public List<string> SequenceCoverageDisplayListWithMods { get; private set; }

		/// <summary>
		/// Visual representation of fragment-level sequence coverage.
		/// Shows coverage from fragment ion identifications (e.g., b/y ions for peptides).
		/// </summary>
		public List<string> FragmentSequenceCoverageDisplayList { get; private set; }

		/// <summary>
		/// Cumulative count of target (non-decoy) groups up to and including this one,
		/// ordered by score. Used for FDR calculation.
		/// </summary>
		public int CumulativeTarget { get; set; }

		/// <summary>
		/// Cumulative count of decoy groups up to and including this one,
		/// ordered by score. Used for FDR calculation.
		/// </summary>
		public int CumulativeDecoy { get; set; }

		/// <summary>
		/// If true, modifications will be displayed in sequence output.
		/// If false, only base sequences will be shown.
		/// </summary>
		public bool DisplayModsOnPeptides { get; set; }

		/// <summary>
		/// Cached string representation of unique sequences for this group.
		/// Populated by GetIdentifiedSequencesOutput method.
		/// </summary>
		private string UniqueSequencesOutput;

		/// <summary>
		/// Cached string representation of sequences shared with other groups.
		/// Populated by GetIdentifiedSequencesOutput method.
		/// </summary>
		private string SharedSequencesOutput;

		#endregion

		#region Methods

		/// <summary>
		/// Prepares output strings for unique and shared sequences identified in this group.
		/// Populates UniqueSequencesOutput and SharedSequencesOutput properties.
		/// </summary>
		/// <param name="labels">Optional SILAC labels for quantification. If provided, sequences will be grouped by label.</param>
		public void GetIdentifiedSequencesOutput(List<SilacLabel> labels)
		{
			var sharedSequences = AllBioPolymersWithSetMods.Except(UniqueBioPolymersWithSetMods);

			if (labels != null)
			{
				// SILAC handling could be implemented here if needed
				// For now, using base implementation
			}

			if (!DisplayModsOnPeptides)
			{
				UniqueSequencesOutput =
					TruncateString(string.Join("|",
						UniqueBioPolymersWithSetMods.Select(p => p.BaseSequence).Distinct()));
				SharedSequencesOutput =
					TruncateString(string.Join("|",
						sharedSequences.Select(p => p.BaseSequence).Distinct()));
			}
			else
			{
				UniqueSequencesOutput =
					TruncateString(string.Join("|",
						UniqueBioPolymersWithSetMods.Select(p => p.FullSequence).Distinct()));
				SharedSequencesOutput =
					TruncateString(string.Join("|",
						sharedSequences.Select(p => p.FullSequence).Distinct()));
			}
		}

		/// <summary>
		/// Returns a tab-separated header line for output files, matching the format of ToString.
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
		/// Format matches the header returned by GetTabSeparatedHeader.
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
			List<double> masses = new List<double>();
			foreach (var sequence in sequences)
			{
				try
				{
					// Calculate mass based on sequence
					// Note: This uses peptide calculation as a fallback; specific biopolymer types
					// may want to override this with their own mass calculation
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

			if (UniqueSequencesOutput != null)
			{
				sb.Append(TruncateString(UniqueSequencesOutput));
			}
			sb.Append("\t");

			if (SharedSequencesOutput != null)
			{
				sb.Append(TruncateString(SharedSequencesOutput));
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

		/// <summary>
		/// Calculates and updates the BioPolymerGroupScore based on PSM scores.
		/// The score is computed as the sum of the best score for each unique base sequence.
		/// </summary>
		public void Score()
		{
			BioPolymerGroupScore = AllPsmsBelowOnePercentFDR.GroupBy(p => p.BaseSequence)
				.Select(p => p.Select(x => x.Score).Max()).Sum();
		}

		/// <summary>
		/// Merges another biopolymer group into this one, combining their members, PSMs, and scores.
		/// Used when groups are determined to represent the same biological entity.
		/// The other group's score is reset to 0 after merging.
		/// </summary>
		/// <param name="otherBioPolymerGroup">The group to merge into this one.</param>
		public void MergeWith(IBioPolymerGroup otherBioPolymerGroup)
		{
			this.BioPolymers.UnionWith(otherBioPolymerGroup.BioPolymers);
			this.AllBioPolymersWithSetMods.UnionWith(otherBioPolymerGroup.AllBioPolymersWithSetMods);
			this.UniqueBioPolymersWithSetMods.UnionWith(otherBioPolymerGroup.UniqueBioPolymersWithSetMods);
			this.AllPsmsBelowOnePercentFDR.UnionWith(otherBioPolymerGroup.AllPsmsBelowOnePercentFDR);
			otherBioPolymerGroup.BioPolymerGroupScore = 0;

			ListOfBioPolymersOrderedByAccession = BioPolymers.OrderBy(p => p.Accession).ToList();
			BioPolymerGroupName = string.Join("|", ListOfBioPolymersOrderedByAccession.Select(p => p.Accession));
		}

		/// <summary>
		/// Creates a new biopolymer group containing only data from a specific file.
		/// The subset group will have the same biopolymers but filtered PSMs, sequences, and intensities.
		/// Used for per-file analysis and output.
		/// </summary>
		/// <param name="fullFilePath">The full path to the file to subset by.</param>
		/// <param name="silacLabels">Optional SILAC labels to apply during subset creation.</param>
		/// <returns>A new BioPolymerGroup containing only data from the specified file.</returns>
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

		#endregion

		#region Equality

		/// <summary>
		/// Determines whether this biopolymer group equals another based on group name.
		/// Two groups are considered equal if they have the same BioPolymerGroupName.
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
		/// Supports comparison with both IBioPolymerGroup and BioPolymerGroup types.
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
		/// Returns a hash code for this biopolymer group based on the group name.
		/// </summary>
		/// <returns>Hash code integer.</returns>
		public override int GetHashCode()
		{
			return BioPolymerGroupName?.GetHashCode() ?? 0;
		}

		#endregion

		#region Private Helpers

		/// <summary>
		/// Truncates a string to MaxStringLength if it exceeds that length.
		/// Returns empty string if input is null or empty.
		/// </summary>
		/// <param name="input">The string to truncate.</param>
		/// <returns>Truncated string or original if within limits.</returns>
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