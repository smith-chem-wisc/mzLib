using MassSpectrometry;
using Omics;
using System.Text;

namespace Proteomics.ProteinGroup
{
    /// <summary>
    /// Represents a group of proteins that share peptides and are grouped together
    /// for quantification and statistical analysis.
    /// </summary>
    /// <typeparam name="TSampleInfo">The sample/file identifier type used to key intensity values.
    /// Must implement <see cref="ISampleInfo"/>.</typeparam>
    public class ProteinGroup<TSampleInfo> where TSampleInfo : ISampleInfo
    {
        /// <summary>
        /// Maximum length for string fields in output. Strings longer than this will be truncated.
        /// Set to 0 or negative to disable truncation.
        /// </summary>
        public static int MaxStringLength { get; set; } = 32000;

        /// <summary>
        /// Creates a new protein group with the specified proteins and peptides.
        /// </summary>
        /// <param name="proteins">The set of proteins that belong to this group.</param>
        /// <param name="peptides">All peptide sequences with set modifications in this group.</param>
        /// <param name="uniquePeptides">Peptides unique to this protein group (not shared with other groups).</param>
        public ProteinGroup(
            HashSet<IBioPolymer> proteins,
            HashSet<IBioPolymerWithSetMods> peptides,
            HashSet<IBioPolymerWithSetMods> uniquePeptides)
        {
            Proteins = proteins;
            ListOfProteinsOrderedByAccession = Proteins.OrderBy(p => p.Accession).ToList();
            ProteinGroupName = string.Join("|", ListOfProteinsOrderedByAccession.Select(p => p.Accession));
            AllPeptides = peptides;
            UniquePeptides = uniquePeptides;
            AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch>();
            SequenceCoverageFraction = new List<double>();
            SequenceCoverageDisplayList = new List<string>();
            SequenceCoverageDisplayListWithMods = new List<string>();
            FragmentSequenceCoverageDisplayList = new List<string>();
            ProteinGroupScore = 0;
            BestPeptideScore = 0;
            BestPeptideQValue = 0;
            QValue = 0;
            IsDecoy = false;
            IsContaminant = false;
            ModsInfo = new List<string>();
            FilesForQuantification = new List<TSampleInfo>();
            IntensitiesByFile = new Dictionary<TSampleInfo, double>();
            CumulativeTarget = 0;
            CumulativeDecoy = 0;
            DisplayModsOnPeptides = false;

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

        /// <summary>
        /// Gets the proteins in this group.
        /// </summary>
        public HashSet<IBioPolymer> Proteins { get; }

        /// <summary>
        /// Gets all proteins in this group ordered alphabetically by accession.
        /// </summary>
        public List<IBioPolymer> ListOfProteinsOrderedByAccession { get; }

        /// <summary>
        /// Gets the protein group name, derived from the accessions of all proteins in the group.
        /// </summary>
        public string ProteinGroupName { get; }

        /// <summary>
        /// Gets all peptides with set modifications in this group.
        /// </summary>
        public HashSet<IBioPolymerWithSetMods> AllPeptides { get; set; }

        /// <summary>
        /// Gets peptides unique to this protein group (not shared with other groups).
        /// </summary>
        public HashSet<IBioPolymerWithSetMods> UniquePeptides { get; set; }

        /// <summary>
        /// Gets all PSMs in this group that pass the 1% FDR threshold.
        /// </summary>
        public HashSet<ISpectralMatch> AllPsmsBelowOnePercentFDR { get; set; }

        /// <summary>
        /// Gets or sets the sequence coverage fraction for each protein in the group.
        /// </summary>
        public List<double> SequenceCoverageFraction { get; set; }

        /// <summary>
        /// Gets or sets the sequence coverage display strings for each protein.
        /// </summary>
        public List<string> SequenceCoverageDisplayList { get; set; }

        /// <summary>
        /// Gets or sets the sequence coverage display strings with modifications for each protein.
        /// </summary>
        public List<string> SequenceCoverageDisplayListWithMods { get; set; }

        /// <summary>
        /// Gets or sets the fragment sequence coverage display strings for each protein.
        /// </summary>
        public List<string> FragmentSequenceCoverageDisplayList { get; set; }

        /// <summary>
        /// Gets or sets the aggregated score for the protein group.
        /// </summary>
        public double ProteinGroupScore { get; set; }

        /// <summary>
        /// Gets or sets the best peptide score in this protein group.
        /// </summary>
        public double BestPeptideScore { get; set; }

        /// <summary>
        /// Gets or sets the best peptide q-value in this protein group.
        /// </summary>
        public double BestPeptideQValue { get; set; }

        /// <summary>
        /// Gets or sets the q-value for this protein group.
        /// </summary>
        public double QValue { get; set; }

        /// <summary>
        /// Gets a value indicating whether this protein group contains decoy proteins.
        /// </summary>
        public bool IsDecoy { get; private set; }

        /// <summary>
        /// Gets a value indicating whether this protein group contains contaminant proteins.
        /// </summary>
        public bool IsContaminant { get; private set; }

        /// <summary>
        /// Gets or sets modification information for this protein group.
        /// </summary>
        public List<string> ModsInfo { get; set; }

        /// <summary>
        /// Gets or sets the list of files/samples used for quantification.
        /// </summary>
        public List<TSampleInfo> FilesForQuantification { get; set; }

        /// <summary>
        /// Gets or sets the dictionary mapping file/sample identifiers to measured intensity values.
        /// </summary>
        public Dictionary<TSampleInfo, double> IntensitiesByFile { get; set; }

        /// <summary>
        /// Gets or sets the cumulative count of target protein groups at this q-value threshold.
        /// Used for FDR calculation.
        /// </summary>
        public int CumulativeTarget { get; set; }

        /// <summary>
        /// Gets or sets the cumulative count of decoy protein groups at this q-value threshold.
        /// Used for FDR calculation.
        /// </summary>
        public int CumulativeDecoy { get; set; }

        /// <summary>
        /// Gets or sets a value indicating whether modifications should be displayed on peptides
        /// in output reports.
        /// </summary>
        public bool DisplayModsOnPeptides { get; set; }

        /// <summary>
        /// Gets or sets the unique peptides output string for display.
        /// </summary>
        public string? UniquePeptidesOutput { get; set; }

        /// <summary>
        /// Gets or sets the shared peptides output string for display.
        /// </summary>
        public string? SharedPeptidesOutput { get; set; }

        /// <summary>
        /// Gets the gene names from all proteins in this group.
        /// </summary>
        public IEnumerable<string> GeneNames => Proteins
            .SelectMany(p => p.GeneNames)
            .Select(g => g.Item2)
            .Distinct();

        /// <summary>
        /// Gets the organisms from all proteins in this group.
        /// </summary>
        public IEnumerable<string> Organisms => Proteins
            .Select(p => p.Organism)
            .Distinct();

        /// <summary>
        /// Gets the number of unique peptides in this protein group.
        /// </summary>
        public int UniquePeptideCount => UniquePeptides.Count;

        /// <summary>
        /// Gets the number of shared peptides in this protein group.
        /// </summary>
        public int SharedPeptideCount => AllPeptides.Count - UniquePeptides.Count;

        /// <summary>
        /// Gets the number of PSMs in this protein group.
        /// </summary>
        public int PsmCount => AllPsmsBelowOnePercentFDR.Count;

        /// <summary>
        /// Truncates a string to <see cref="MaxStringLength"/> if it exceeds that length.
        /// Returns the original string if MaxStringLength is 0 or negative (disabled).
        /// </summary>
        /// <param name="input">The string to truncate.</param>
        /// <returns>The truncated string, or the original if truncation is disabled or not needed.</returns>
        private static string TruncateString(string? input)
        {
            if (string.IsNullOrEmpty(input))
                return string.Empty;

            if (MaxStringLength <= 0 || input.Length <= MaxStringLength)
                return input;

            return input.Substring(0, MaxStringLength);
        }

        /// <summary>
        /// Returns a tab-separated header line for output files.
        /// </summary>
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

            if (FilesForQuantification != null && FilesForQuantification.Count > 0)
            {
                bool unfractionated = FilesForQuantification.Select(p => p.Fraction).Distinct().Count() == 1;
                bool conditionsUndefined = FilesForQuantification.All(p => string.IsNullOrEmpty(p.Condition));

                // Check if this is SILAC-labeled data by looking for non-existent file paths
                bool silacExperimentalDesign =
                    FilesForQuantification.Any(p => !string.IsNullOrEmpty(p.FullFilePathWithExtension)
                        && !File.Exists(p.FullFilePathWithExtension));

                foreach (var sampleGroup in FilesForQuantification.GroupBy(p => p.Condition))
                {
                    foreach (var sample in sampleGroup.GroupBy(p => p.BiologicalReplicate).OrderBy(p => p.Key))
                    {
                        var firstSample = sample.First();

                        if ((conditionsUndefined && unfractionated) || silacExperimentalDesign)
                        {
                            // Use the display name as the intensity header
                            sb.Append("Intensity_" + firstSample.DisplayName + "\t");
                        }
                        else
                        {
                            // Label the header with condition and biorep number
                            sb.Append("Intensity_" + firstSample.Condition + "_" +
                                      (firstSample.BiologicalReplicate + 1) + "\t");
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

        /// <summary>
        /// Returns a tab-separated string representation of this protein group for output files.
        /// </summary>
        public override string ToString()
        {
            var sb = new StringBuilder();

            // list of protein accession numbers
            sb.Append(ProteinGroupName);
            sb.Append("\t");

            // genes
            sb.Append(TruncateString(string.Join("|",
                ListOfProteinsOrderedByAccession.Select(p => p.GeneNames.Select(x => x.Item2).FirstOrDefault()))));
            sb.Append("\t");

            // organisms
            sb.Append(TruncateString(string.Join("|",
                ListOfProteinsOrderedByAccession.Select(p => p.Organism).Distinct())));
            sb.Append("\t");

            // list of protein names
            sb.Append(TruncateString(string.Join("|",
                ListOfProteinsOrderedByAccession.Select(p => p.FullName).Distinct())));
            sb.Append("\t");

            // list of masses
            var sequences = ListOfProteinsOrderedByAccession.Select(p => p.BaseSequence).Distinct();
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

            // number of proteins in group
            sb.Append("" + Proteins.Count);
            sb.Append("\t");

            // list of unique peptides
            if (UniquePeptidesOutput != null)
            {
                sb.Append(TruncateString(UniquePeptidesOutput));
            }

            sb.Append("\t");

            // list of shared peptides
            if (SharedPeptidesOutput != null)
            {
                sb.Append(TruncateString(SharedPeptidesOutput));
            }

            sb.Append("\t");

            // number of peptides
            if (!DisplayModsOnPeptides)
            {
                sb.Append("" + AllPeptides.Select(p => p.BaseSequence).Distinct().Count());
            }
            else
            {
                sb.Append("" + AllPeptides.Select(p => p.FullSequence).Distinct().Count());
            }

            sb.Append("\t");

            // number of unique peptides
            if (!DisplayModsOnPeptides)
            {
                sb.Append("" + UniquePeptides.Select(p => p.BaseSequence).Distinct().Count());
            }
            else
            {
                sb.Append("" + UniquePeptides.Select(p => p.FullSequence).Distinct().Count());
            }

            sb.Append("\t");

            // sequence coverage percent
            sb.Append(TruncateString(string.Join("|",
                SequenceCoverageFraction.Select(p => string.Format("{0:0.#####}", p)))));
            sb.Append("\t");

            // sequence coverage
            sb.Append(TruncateString(string.Join("|", SequenceCoverageDisplayList)));
            sb.Append("\t");

            // sequence coverage with mods
            sb.Append(TruncateString(string.Join("|", SequenceCoverageDisplayListWithMods)));
            sb.Append("\t");

            // fragment sequence coverage
            sb.Append(TruncateString(string.Join("|", FragmentSequenceCoverageDisplayList)));
            sb.Append("\t");

            // Detailed mods information list
            sb.Append(TruncateString(string.Join("|", ModsInfo)));
            sb.Append("\t");

            // MS1 intensity (retrieved from FlashLFQ in the SearchTask)
            if (IntensitiesByFile != null && FilesForQuantification != null)
            {
                foreach (var sampleGroup in FilesForQuantification.GroupBy(p => p.Condition))
                {
                    foreach (var sample in sampleGroup.GroupBy(p => p.BiologicalReplicate).OrderBy(p => p.Key))
                    {
                        // if the samples are fractionated, the protein will only have 1 intensity in the first fraction
                        // and the other fractions will be zero. we could find the first/only fraction with an intensity,
                        // but simply summing the fractions is easier than finding the single non-zero value
                        double summedIntensity = sample.Sum(file => IntensitiesByFile.TryGetValue(file, out var intensity) ? intensity : 0);

                        if (summedIntensity > 0)
                        {
                            sb.Append(summedIntensity);
                        }

                        sb.Append("\t");
                    }
                }
            }

            // number of PSMs for listed peptides
            sb.Append("" + AllPsmsBelowOnePercentFDR.Count);
            sb.Append("\t");

            // isDecoy
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

            // cumulative target
            sb.Append(CumulativeTarget);
            sb.Append("\t");

            // cumulative decoy
            sb.Append(CumulativeDecoy);
            sb.Append("\t");

            // q value
            sb.Append(QValue);
            sb.Append("\t");

            // best peptide score
            sb.Append(BestPeptideScore);
            sb.Append("\t");

            // best peptide q value
            sb.Append(BestPeptideQValue);

            return sb.ToString();
        }
        // this method is only used internally, to make protein grouping faster
        // this is NOT an output and is NOT used for protein FDR calculations
        public void Score()
        {
            // sum the scores of the best PSM per base sequence
            ProteinGroupScore = AllPsmsBelowOnePercentFDR.GroupBy(p => p.BaseSequence)
                .Select(p => p.Select(x => x.Score).Max()).Sum();
        }
        public void CalculateSequenceCoverage()
        {
            var proteinsWithUnambigSeqPsms = new Dictionary<IBioPolymer, List<IBioPolymerWithSetMods>>();
            var proteinsWithPsmsWithLocalizedMods = new Dictionary<IBioPolymer, List<IBioPolymerWithSetMods>>();

            foreach (var protein in Proteins)
            {
                proteinsWithUnambigSeqPsms.Add(protein, new List<IBioPolymerWithSetMods>());
                proteinsWithPsmsWithLocalizedMods.Add(protein, new List<IBioPolymerWithSetMods>());
            }

            foreach (var psm in AllPsmsBelowOnePercentFDR)
            {
                // null BaseSequence means that the amino acid sequence is ambiguous; do not use these to calculate sequence coverage
                if (psm.BaseSequence != null)
                {
                    psm.GetAminoAcidCoverage();

                    foreach (var peptide in psm.BestMatchingBioPolymersWithSetMods.Select(psm => psm.SpecificBioPolymer).DistinctBy(pep => pep.FullSequence))
                    {
                        // might be unambiguous but also shared; make sure this protein group contains this peptide+protein combo
                        if (Proteins.Contains(peptide.Parent))
                        {
                            proteinsWithUnambigSeqPsms[peptide.Parent].Add(peptide);

                            // null FullSequence means that mods were not successfully localized; do not display them on the sequence coverage mods info
                            if (peptide.FullSequence != null)
                            {
                                proteinsWithPsmsWithLocalizedMods[peptide.Parent].Add(peptide);
                            }
                        }
                    }

                }
            }

            //Calculate sequence coverage at the amino acid level by looking at fragment specific coverage
            //loop through proteins
            foreach (IBioPolymer protein in ListOfProteinsOrderedByAccession)
            {
                //create a hash set for storing covered one-based residue numbers of protein
                HashSet<int> coveredResiduesInProteinOneBased = new();

                //loop through PSMs
                foreach (SpectralMatch psm in AllPsmsBelowOnePercentFDR.Where(psm => psm.BaseSequence != null))
                {
                    //Calculate the covered bases within the psm. This is one based numbering for the peptide only
                    psm.GetAminoAcidCoverage();
                    if (psm.FragmentCoveragePositionInPeptide == null) continue;
                    //loop through each peptide within the psm
                    IEnumerable<IBioPolymerWithSetMods> pwsms = psm.BestMatchingBioPolymersWithSetMods.Select(p => p.SpecificBioPolymer)
                        .Where(p => p.Parent.Accession == protein.Accession);
                    foreach (var pwsm in pwsms)
                    {
                        //create a hashset to store the covered residues for the peptide, converted to the corresponding indices of the protein
                        HashSet<int> coveredResiduesInPeptide = new();
                        //add the peptide start position within the protein to each covered index of the psm
                        foreach (var position in psm.FragmentCoveragePositionInPeptide)
                        {
                            coveredResiduesInPeptide.Add(position + pwsm.OneBasedStartResidue -
                                                         1); //subtract one because these are both one based
                        }

                        //Add the peptide specific positions, to the overall hashset for the protein
                        coveredResiduesInProteinOneBased.UnionWith(coveredResiduesInPeptide);
                    }
                }

                // create upper/lowercase string
                char[] fragmentCoverageArray = protein.BaseSequence.ToLower().ToCharArray();
                foreach (var residue in coveredResiduesInProteinOneBased)
                {
                    fragmentCoverageArray[residue - 1] = char.ToUpper(fragmentCoverageArray[residue - 1]);
                }

                FragmentSequenceCoverageDisplayList.Add(new string(fragmentCoverageArray));
            }

            //Calculates the coverage at the peptide level... if a peptide is present all of the AAs in the peptide are covered
            foreach (var protein in ListOfProteinsOrderedByAccession)
            {
                HashSet<int> coveredOneBasedResidues = new HashSet<int>();

                // get residue numbers of each peptide in the protein and identify them as observed if the sequence is unambiguous
                foreach (var peptide in proteinsWithUnambigSeqPsms[protein])
                {
                    for (int i = peptide.OneBasedStartResidue; i <= peptide.OneBasedEndResidue; i++)
                    {
                        coveredOneBasedResidues.Add(i);
                    }
                }

                // calculate sequence coverage percent
                double seqCoverageFract = (double)coveredOneBasedResidues.Count / protein.Length;

                // add the percent coverage
                SequenceCoverageFraction.Add(seqCoverageFract);

                // convert the observed amino acids to upper case if they are unambiguously observed
                string sequenceCoverageDisplay = protein.BaseSequence.ToLower();
                var coverageArray = sequenceCoverageDisplay.ToCharArray();
                foreach (var obsResidueLocation in coveredOneBasedResidues)
                {
                    coverageArray[obsResidueLocation - 1] = char.ToUpper(coverageArray[obsResidueLocation - 1]);
                }

                sequenceCoverageDisplay = new string(coverageArray);

                // add the coverage display
                SequenceCoverageDisplayList.Add(sequenceCoverageDisplay);

                // put mods in the sequence coverage display
                // get mods to display in sequence (only unambiguously identified mods)
                var modsOnThisProtein = new HashSet<KeyValuePair<int, Modification>>();
                foreach (var pep in proteinsWithPsmsWithLocalizedMods[protein])
                {
                    foreach (var mod in pep.AllModsOneIsNterminus)
                    {
                        if (!mod.Value.ModificationType.Contains("PeptideTermMod")
                            && !mod.Value.ModificationType.Contains("Common Variable")
                            && !mod.Value.ModificationType.Contains("Common Fixed"))
                        {
                            modsOnThisProtein.Add(
                                new KeyValuePair<int, Modification>(pep.OneBasedStartResidue + mod.Key - 2,
                                    mod.Value));
                        }
                    }
                }

                var tempMods = modsOnThisProtein.OrderBy(p => p.Key).ToList();
                foreach (var mod in tempMods)
                {
                    if (mod.Value.LocationRestriction.Equals("N-terminal."))
                    {
                        sequenceCoverageDisplay = sequenceCoverageDisplay.Insert(
                            0,
                            $"[{mod.Value.IdWithMotif}]-");
                    }
                    else if (mod.Value.LocationRestriction.Equals("Anywhere."))
                    {
                        int modStringIndex = sequenceCoverageDisplay.Length - (protein.Length - mod.Key);
                        sequenceCoverageDisplay = sequenceCoverageDisplay.Insert(
                            modStringIndex,
                            $"[{mod.Value.IdWithMotif}]");
                    }
                    else if (mod.Value.LocationRestriction.Equals("C-terminal."))
                    {
                        sequenceCoverageDisplay = sequenceCoverageDisplay.Insert(
                            sequenceCoverageDisplay.Length,
                            $"-[{mod.Value.IdWithMotif}]");
                    }
                }

                SequenceCoverageDisplayListWithMods.Add(sequenceCoverageDisplay);

                if (!modsOnThisProtein.Any())
                {
                    continue;
                }

                // calculate spectral count % of modified observations
                var pepModTotals = new List<int>(); // count of modified peptides for each mod/index
                var pepTotals = new List<int>(); // count of all peptides for each mod/index
                var modIndex = new List<(int index, string modName)>(); // index and name of the modified position

                foreach (var pep in proteinsWithPsmsWithLocalizedMods[protein])
                {
                    foreach (var mod in pep.AllModsOneIsNterminus)
                    {
                        int pepNumTotal = 0; //For one mod, The total Pep Num

                        if (mod.Value.ModificationType.Contains("Common Variable")
                            || mod.Value.ModificationType.Contains("Common Fixed")
                            || mod.Value.LocationRestriction.Equals(ModLocationOnPeptideOrProtein.PepC)
                            || mod.Value.LocationRestriction.Equals(ModLocationOnPeptideOrProtein.NPep))
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
                            indexInProtein = pep.OneBasedStartResidue + mod.Key - 2;
                        }
                        else if (mod.Value.LocationRestriction.Equals("C-terminal."))
                        {
                            indexInProtein = protein.Length;
                        }
                        else
                        {
                            // In case it's a peptide terminal mod, skip!
                            // we don't want this annotated in the protein's modifications
                            continue;
                        }

                        var modKey = (indexInProtein, mod.Value.IdWithMotif);
                        if (modIndex.Contains(modKey))
                        {
                            pepModTotals[modIndex.IndexOf(modKey)] += 1;
                        }
                        else
                        {
                            modIndex.Add(modKey);
                            foreach (var pept in proteinsWithPsmsWithLocalizedMods[protein])
                            {
                                if (indexInProtein >= pept.OneBasedStartResidue - (indexInProtein == 1 ? 1 : 0)
                                    && indexInProtein <= pept.OneBasedEndResidue)
                                {
                                    pepNumTotal += 1;
                                }
                            }

                            pepTotals.Add(pepNumTotal);
                            pepModTotals.Add(1);
                        }
                    }
                }

                var modStrings = new List<(int aaNum, string part)>();
                for (int i = 0; i < pepModTotals.Count; i++)
                {
                    string aa = modIndex[i].index.ToString();
                    string modName = modIndex[i].modName.ToString();
                    string occupancy = ((double)pepModTotals[i] / (double)pepTotals[i]).ToString("F2");
                    string fractOccupancy = $"{pepModTotals[i].ToString()}/{pepTotals[i].ToString()}";
                    string tempString = ($"#aa{aa}[{modName},info:occupancy={occupancy}({fractOccupancy})]");
                    modStrings.Add((modIndex[i].index, tempString));
                }

                var modInfoString = string.Join(";", modStrings.OrderBy(x => x.aaNum).Select(x => x.part));

                if (!string.IsNullOrEmpty(modInfoString))
                {
                    ModsInfo.Add(modInfoString);
                }
            }
        }
        public void MergeProteinGroupWith(ProteinGroup other)
        {
            this.Proteins.UnionWith(other.Proteins);
            this.AllPeptides.UnionWith(other.AllPeptides);
            this.UniquePeptides.UnionWith(other.UniquePeptides);
            this.AllPsmsBelowOnePercentFDR.UnionWith(other.AllPsmsBelowOnePercentFDR);
            other.ProteinGroupScore = 0;

            ListOfProteinsOrderedByAccession = Proteins.OrderBy(p => p.Accession).ToList();

            ProteinGroupName = string.Join("|", ListOfProteinsOrderedByAccession.Select(p => p.Accession));
        }
        public ProteinGroup ConstructSubsetProteinGroup(string fullFilePath, List<SilacLabel> silacLabels = null)
        {
            var allPsmsForThisFile =
                new HashSet<SpectralMatch>(
                    AllPsmsBelowOnePercentFDR.Where(p => p.FullFilePath.Equals(fullFilePath)));
            var allPeptidesForThisFile =
                new HashSet<IBioPolymerWithSetMods>(
                    allPsmsForThisFile.SelectMany(p => p.BestMatchingBioPolymersWithSetMods.Select(v => v.SpecificBioPolymer)));
            var allUniquePeptidesForThisFile =
                new HashSet<IBioPolymerWithSetMods>(UniquePeptides.Intersect(allPeptidesForThisFile));

            ProteinGroup subsetPg = new ProteinGroup(Proteins, allPeptidesForThisFile, allUniquePeptidesForThisFile)
            {
                AllPsmsBelowOnePercentFDR = allPsmsForThisFile,
                DisplayModsOnPeptides = DisplayModsOnPeptides
            };

            SpectraFileInfo spectraFileInfo = null;
            if (FilesForQuantification != null)
            {
                spectraFileInfo = FilesForQuantification.Where(p => p.FullFilePathWithExtension == fullFilePath)
                    .FirstOrDefault();
                //check that file name wasn't changed (can occur in SILAC searches)
                if (!silacLabels.IsNullOrEmpty() && spectraFileInfo == null)
                {
                    foreach (SilacLabel label in silacLabels)
                    {
                        string fakeFilePath = SilacConversions
                            .GetHeavyFileInfo(new SpectraFileInfo(fullFilePath, "", 0, 0, 0), label)
                            .FullFilePathWithExtension;
                        spectraFileInfo = FilesForQuantification.Where(p => p.FullFilePathWithExtension == fakeFilePath)
                            .FirstOrDefault();
                        if (spectraFileInfo != null)
                        {
                            break;
                        }
                    }

                    //if still no hits, might be SILAC turnover
                    if (spectraFileInfo == null)
                    {
                        string filepathWithoutExtension = Path.Combine(Path.GetDirectoryName(fullFilePath),
                            Path.GetFileNameWithoutExtension(fullFilePath));
                        string extension = Path.GetExtension(fullFilePath);
                        string fakeFilePath = filepathWithoutExtension + SilacConversions.ORIGINAL_TURNOVER_LABEL_NAME +
                                              extension;
                        spectraFileInfo = FilesForQuantification.Where(p => p.FullFilePathWithExtension == fakeFilePath)
                            .FirstOrDefault();
                    }
                }

                subsetPg.FilesForQuantification = new List<SpectraFileInfo> { spectraFileInfo };
            }

            if (IntensitiesByFile == null)
            {
                subsetPg.IntensitiesByFile = null;
            }
            else
            {
                subsetPg.IntensitiesByFile = new Dictionary<SpectraFileInfo, double>
                    { { spectraFileInfo, IntensitiesByFile[spectraFileInfo] } };
            }

            return subsetPg;
        }
        public bool Equals(ProteinGroup grp)
        {
            //Check for null and compare run-time types.
            if (grp == null)
            {
                return false;
            }
            else if (!this.ListOfProteinsOrderedByAccession.Select(a => a.Accession).ToList().SequenceEqual(grp.ListOfProteinsOrderedByAccession.Select(a => a.Accession).ToList()))
            {
                return false;
            }

            return true;
        }
    }

    /// <summary>
    /// Protein group for label-free quantification using <see cref="SpectraFileInfo"/> as the file key.
    /// </summary>
    public class LabelFreeProteinGroup : ProteinGroup<SpectraFileInfo>
    {
        /// <summary>
        /// Creates a new label-free protein group.
        /// </summary>
        public LabelFreeProteinGroup(
            HashSet<IBioPolymer> proteins,
            HashSet<IBioPolymerWithSetMods> peptides,
            HashSet<IBioPolymerWithSetMods> uniquePeptides)
            : base(proteins, peptides, uniquePeptides)
        {
        }
    }

    /// <summary>
    /// Protein group for isobaric (TMT/iTRAQ) quantification using <see cref="IIsobaricQuantSampleInfo"/> as the sample key.
    /// </summary>
    public class IsobaricProteinGroup : ProteinGroup<IIsobaricQuantSampleInfo>
    {
        /// <summary>
        /// Creates a new isobaric protein group.
        /// </summary>
        public IsobaricProteinGroup(
            HashSet<IBioPolymer> proteins,
            HashSet<IBioPolymerWithSetMods> peptides,
            HashSet<IBioPolymerWithSetMods> uniquePeptides)
            : base(proteins, peptides, uniquePeptides)
        {
        }
    }
}