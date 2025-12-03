using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using Omics.Modifications;
using System.Text;

namespace Omics.BioPolymerGroup
{
    public class BioPolymerGroup : IBioPolymerGroup
    {
        // Flags
        public bool IsDecoy { get; init; }
        public bool IsContaminant { get; init; }

        // Core collections
        public HashSet<IBioPolymer> BioPolymers { get; set; } = new();
        public string BioPolymerGroupName { get; private set; } = string.Empty;
        public double BioPolymerGroupScore { get; private set; }
        public List<SpectraFileInfo> FilesForQuantification { get; set; }

        // BioPolymerWithSetMods
        public HashSet<IBioPolymerWithSetMods> AllBioPolymerWithSetMods { get; set; } = new();
        public HashSet<IBioPolymerWithSetMods> UniqueBioPolymerWithSetMods { get; set; } = new();

        // Coverage and display
        public List<double> SequenceCoverageFraction { get; } = new();
        public List<string> SequenceCoverageDisplayList { get; } = new();
        public List<string> SequenceCoverageDisplayListWithMods { get; } = new();
        public List<string> FragmentSequenceCoverageDisplayList { get; } = new();

        // FDR/Q-value tracking
        public double QValue { get; set; }
        public double BestBioPolymerWithSetModQValue { get; set; }
        public double BestBioPolymerWithSetModScore { get; set; }
        public int CumulativeTarget { get; set; }
        public int CumulativeDecoy { get; set; }
        public bool DisplayModsOnBioPolymerWithSetMods { get; set; }

        // Mods aggregate info
        public List<string> ModsInfo { get; } = new();

        // Convenience views
        public List<IBioPolymer> ListOfBioPolymersOrderedByAccession { get; private set; } =
            BioPolymers.OrderBy(b => b.Accession, StringComparer.Ordinal).ToList();

        public string UniqueBioPolymerWithSetModsOutput =>
            FormatBioPolymerWithSetModOutput(UniqueBioPolymerWithSetMods);

        public string SharedBioPolymerWithSetModsOutput =>
            FormatBioPolymerWithSetModOutput(AllBioPolymerWithSetMods.Except(UniqueBioPolymerWithSetMods).ToHashSet());

        // Ctors
        public BioPolymerGroup() { }

        public BioPolymerGroup(
            string groupName,
            IEnumerable<IBioPolymer>? biopolymers = null,
            bool isDecoy = false,
            bool isContaminant = false)
        {
            BioPolymerGroupName = groupName ?? string.Empty;
            if (biopolymers != null) BioPolymers = new HashSet<IBioPolymer>(biopolymers);
            IsDecoy = isDecoy;
            IsContaminant = isContaminant;
        }
        // neccessary External Values
        public int MaxLengthOfOutput { get; set; } = int.MaxValue;

        // Equality by group identity (name + decoy/contaminant flags)
        private static string CanonicalBioPolymerKey(HashSet<IBioPolymer> biopolymers) =>
            string.Join("|", biopolymers.Select(b => b.Accession).OrderBy(a => a, StringComparer.Ordinal));

        public override bool Equals(object? obj) => Equals(obj as IBioPolymerGroup);

        public override int GetHashCode()
        {
            return HashCode.Combine(
                StringComparer.Ordinal.GetHashCode(BioPolymerGroupName ?? string.Empty),
                IsDecoy,
                IsContaminant,
                StringComparer.Ordinal.GetHashCode(CanonicalBioPolymerKey(BioPolymers)));
        }

        // Helpers
        private string FormatBioPolymerWithSetModOutput(HashSet<IBioPolymerWithSetMods> peptides)
        {
            if (peptides == null || peptides.Count == 0) return string.Empty;

            // Output as semicolon-separated base sequences; include mods if requested
            var parts = peptides
                .OrderBy(p => p.BaseSequence, StringComparer.Ordinal)
                .Select(p => DisplayModsOnBioPolymerWithSetMods ? p.FullSequence : p.BaseSequence);

            return string.Join(";", parts);
        }
        public string GetTabSeparatedHeader()
        {
            var sb = new StringBuilder();
            sb.Append("BioPolymer Accession" + '\t');
            sb.Append("Gene" + '\t');
            sb.Append("Organism" + '\t');
            sb.Append("BioPolymer Full Name" + '\t');
            sb.Append("BioPolymer Unmodified Mass" + '\t');
            sb.Append("Number of BioPolymers in Group" + '\t');
            sb.Append("Unique BioPolymerWithSetMods" + '\t');
            sb.Append("Shared BioPolymerWithSetMods" + '\t');
            sb.Append("Number of BioPolymerWithSetMods" + '\t');
            sb.Append("Number of Unique BioPolymerWithSetMods" + '\t');
            sb.Append("Sequence Coverage Fraction" + '\t');
            sb.Append("Sequence Coverage" + '\t');
            sb.Append("Sequence Coverage with Mods" + '\t');
            sb.Append("Fragment Sequence Coverage" + '\t');
            sb.Append("Modification Info List" + "\t");
            if (FilesForQuantification != null)
            {
                bool unfractionated = FilesForQuantification.Select(p => p.Fraction).Distinct().Count() == 1;
                bool conditionsUndefined = FilesForQuantification.All(p => string.IsNullOrEmpty(p.Condition));

                // this is a hacky way to test for SILAC-labeled data...
                // Currently SILAC will report 1 column of intensities per label per spectra file, and is NOT summarized
                // into biorep-level intensity values. the SILAC code uses the "condition" field to organize this info,
                // even if the experimental design is not defined by the user. So the following bool is a way to distinguish
                // between experimental design being used in SILAC automatically vs. being defined by the user
                bool silacExperimentalDesign =
                    FilesForQuantification.Any(p => !File.Exists(p.FullFilePathWithExtension));

                foreach (var sampleGroup in FilesForQuantification.GroupBy(p => p.Condition))
                {
                    foreach (var sample in sampleGroup.GroupBy(p => p.BiologicalReplicate).OrderBy(p => p.Key))
                    {
                        if ((conditionsUndefined && unfractionated) || silacExperimentalDesign)
                        {
                            // if the data is unfractionated and the conditions haven't been defined, just use the file name as the intensity header
                            sb.Append("Intensity_" + sample.First().FilenameWithoutExtension + "\t");
                        }
                        else
                        {
                            // if the data is fractionated and/or the conditions have been defined, label the header w/ the condition and biorep number
                            sb.Append("Intensity_" + sample.First().Condition + "_" +
                                      (sample.First().BiologicalReplicate + 1) + "\t");
                        }
                    }
                }
            }

            sb.Append("Number of Spectrum Matches" + '\t');
            sb.Append("BioPolymer Decoy/Contaminant/Target" + '\t');
            sb.Append("BioPolymer Cumulative Target" + '\t');
            sb.Append("BioPolymer Cumulative Decoy" + '\t');
            sb.Append("BioPolymer QValue" + '\t');
            sb.Append("Best BioPolymerWithSetMods Score" + '\t');
            sb.Append("Best BioPolymerWithSetMods Notch QValue");
            return sb.ToString();
        }
        public override string ToString()
        {
            var sb = new StringBuilder();

            // list of protein accession numbers
            sb.Append(BioPolymerGroupName);
            sb.Append("\t");

            // genes
            sb.Append(string.Join("|",
                ListOfBioPolymersOrderedByAccession.Select(p => p.GeneNames.Select(x => x.Item2).FirstOrDefault())).Substring(0,MaxLengthOfOutput));
            sb.Append("\t");

            // organisms
            sb.Append(string.Join("|",
                ListOfBioPolymersOrderedByAccession.Select(p => p.Organism).Distinct()).Substring(0, MaxLengthOfOutput));
            sb.Append("\t");

            // list of protein names
            sb.Append(string.Join("|",
                ListOfBioPolymersOrderedByAccession.Select(p => p.FullName).Distinct()).Substring(0, MaxLengthOfOutput));
            sb.Append("\t");

            // list of masses
            var sequences = ListOfBioPolymersOrderedByAccession.Select(p => p.BaseSequence).Distinct();
            List<double> masses = new List<double>();
            foreach (var sequence in sequences)
            {
                try
                {
                    if (GlobalVariables.AnalyteType == AnalyteType.Oligo)
                        masses.Add(new OligoWithSetMods(sequence, GlobalVariables.AllRnaModsKnownDictionary).MonoisotopicMass);
                    else
                        masses.Add(new Proteomics.AminoAcidPolymer.BioPolymerWithSetMod(sequence).MonoisotopicMass);

                }
                catch (System.Exception)
                {
                    masses.Add(double.NaN);
                }
            }

            sb.Append(string.Join("|", masses).Substring(0, MaxLengthOfOutput));
            sb.Append("\t");

            // number of proteins in group
            sb.Append("" + BioPolymers.Count);
            sb.Append("\t");

            // list of unique peptides
            if (UniqueBioPolymerWithSetModsOutput != null)
            {
                sb.Append(UniqueBioPolymerWithSetModsOutput.Substring(0, MaxLengthOfOutput));
            }

            sb.Append("\t");

            // list of shared biopolymers
            if (SharedBioPolymerWithSetModsOutput != null)
            {
                sb.Append(SharedBioPolymerWithSetModsOutput.Substring(0, MaxLengthOfOutput));
            }

            sb.Append("\t");

            // number of biopolymers
            if (!DisplayModsOnBioPolymerWithSetMods)
            {
                sb.Append("" + AllBioPolymerWithSetMods.Select(p => p.BaseSequence).Distinct().Count());
            }
            else
            {
                sb.Append("" + AllBioPolymerWithSetMods.Select(p => p.FullSequence).Distinct().Count());
            }

            sb.Append("\t");

            // number of unique biopolymers
            if (!DisplayModsOnBioPolymerWithSetMods)
            {
                sb.Append("" + UniqueBioPolymerWithSetMods.Select(p => p.BaseSequence).Distinct().Count());
            }
            else
            {
                sb.Append("" + UniqueBioPolymerWithSetMods.Select(p => p.FullSequence).Distinct().Count());
            }

            sb.Append("\t");

            // sequence coverage percent
            sb.Append(string.Join("|",
                SequenceCoverageFraction.Select(p => string.Format("{0:0.#####}", p))).Substring(0, MaxLengthOfOutput));
            sb.Append("\t");

            // sequence coverage
            sb.Append(string.Join("|", SequenceCoverageDisplayList).Substring(0, MaxLengthOfOutput));
            sb.Append("\t");

            // sequence coverage with mods
            sb.Append(string.Join("|", SequenceCoverageDisplayListWithMods).Substring(0, MaxLengthOfOutput));
            sb.Append("\t");

            // fragment sequence coverage
            sb.Append(string.Join("|", FragmentSequenceCoverageDisplayList).Substring(0, MaxLengthOfOutput));
            sb.Append("\t");

            //Detailed mods information list
            sb.Append(string.Join("|", ModsInfo).Substring(0, MaxLengthOfOutput));
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
                        double summedIntensity = sample.Sum(file => IntensitiesByFile[file]);

                        if (summedIntensity > 0)
                        {
                            sb.Append(summedIntensity);
                        }

                        sb.Append("\t");
                    }
                }
            }

            // number of PSMs for listed biopolymers
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

            // best biopolymer score
            sb.Append(BestBioPolymerWithSetModScore);
            sb.Append("\t");

            // best biopolymer q value
            sb.Append(BestBioPolymerWithSetModQValue);

            return sb.ToString();
        }
        // this method is only used internally, to make biopolymer grouping faster
        // this is NOT an output and is NOT used for biopolymer FDR calculations
        public void Score()
        {
            // sum the scores of the best PSM per base sequence
            BioPolymerGroupScore = AllPsmsBelowOnePercentFDR.GroupBy(p => p.BaseSequence)
                .Select(p => p.Select(x => x.Score).Max()).Sum();
        }
        public void CalculateSequenceCoverage()
        {
            var proteinsWithUnambigSeqPsms = new Dictionary<IBioPolymer, List<IBioPolymerWithSetMods>>();
            var proteinsWithPsmsWithLocalizedMods = new Dictionary<IBioPolymer, List<IBioPolymerWithSetMods>>();

            foreach (var biopolymer in BioPolymers)
            {
                proteinsWithUnambigSeqPsms.Add(biopolymer, new List<IBioPolymerWithSetMods>());
                proteinsWithPsmsWithLocalizedMods.Add(biopolymer, new List<IBioPolymerWithSetMods>());
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
                        if (BioPolymers.Contains(peptide.Parent))
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
            foreach (IBioPolymer protein in ListOfBioPolymersOrderedByAccession)
            {
                //create a hash set for storing covered one-based residue numbers of protein
                HashSet<int> coveredResiduesInProteinOneBased = new();

                //loop through PSMs
                foreach (SpectralMatch psm in AllPsmsBelowOnePercentFDR.Where(psm => psm.BaseSequence != null))
                {
                    //Calculate the covered bases within the psm. This is one based numbering for the peptide only
                    psm.GetAminoAcidCoverage();
                    if (psm.FragmentCoveragePositionInBioPolymerWithSetMod == null) continue;
                    //loop through each peptide within the psm
                    IEnumerable<IBioPolymerWithSetMods> pwsms = psm.BestMatchingBioPolymersWithSetMods.Select(p => p.SpecificBioPolymer)
                        .Where(p => p.Parent.Accession == protein.Accession);
                    foreach (var pwsm in pwsms)
                    {
                        //create a hashset to store the covered residues for the peptide, converted to the corresponding indices of the protein
                        HashSet<int> coveredResiduesInBioPolymerWithSetMod = new();
                        //add the peptide start position within the protein to each covered index of the psm
                        foreach (var position in psm.FragmentCoveragePositionInBioPolymerWithSetMod)
                        {
                            coveredResiduesInBioPolymerWithSetMod.Add(position + pwsm.OneBasedStartResidue -
                                                         1); //subtract one because these are both one based
                        }

                        //Add the peptide specific positions, to the overall hashset for the protein
                        coveredResiduesInProteinOneBased.UnionWith(coveredResiduesInBioPolymerWithSetMod);
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
            foreach (var protein in ListOfBioPolymersOrderedByAccession)
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
                        if (!mod.Value.ModificationType.Contains("BioPolymerWithSetModTermMod")
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
                            || mod.Value.LocationRestriction.Equals(ModLocationOnBioPolymerWithSetModOrProtein.PepC)
                            || mod.Value.LocationRestriction.Equals(ModLocationOnBioPolymerWithSetModOrProtein.NPep))
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
        public void MergeWith(BioPolymerGroup other)
        {
            this.BioPolymers.UnionWith(other.BioPolymers);
            this.AllBioPolymerWithSetMods.UnionWith(other.AllBioPolymerWithSetMods);
            this.UniqueBioPolymerWithSetMods.UnionWith(other.UniqueBioPolymerWithSetMods);
            this.AllPsmsBelowOnePercentFDR.UnionWith(other.AllPsmsBelowOnePercentFDR);
            other.BioPolymerGroupScore = 0;

            ListOfBioPolymersOrderedByAccession = BioPolymers.OrderBy(p => p.Accession).ToList();

            BioPolymerGroupName = string.Join("|", ListOfBioPolymersOrderedByAccession.Select(p => p.Accession));
        }
        public BioPolymerGroup ConstructSubsetBioPolymerGroup(string fullFilePath, List<SilacLabel> silacLabels = null)
        {
            var allPsmsForThisFile =
                new HashSet<SpectralMatch>(
                    AllPsmsBelowOnePercentFDR.Where(p => p.FullFilePath.Equals(fullFilePath)));
            var allBioPolymerWithSetModsForThisFile =
                new HashSet<IBioPolymerWithSetMods>(
                    allPsmsForThisFile.SelectMany(p => p.BestMatchingBioPolymersWithSetMods.Select(v => v.SpecificBioPolymer)));
            var allUniqueBioPolymerWithSetModsForThisFile =
                new HashSet<IBioPolymerWithSetMods>(UniqueBioPolymerWithSetMods.Intersect(allBioPolymerWithSetModsForThisFile));

            BioPolymerGroup subsetPg = new BioPolymerGroup(BioPolymers, allBioPolymerWithSetModsForThisFile, allUniqueBioPolymerWithSetModsForThisFile)
            {
                AllPsmsBelowOnePercentFDR = allPsmsForThisFile,
                DisplayModsOnBioPolymerWithSetMods = DisplayModsOnBioPolymerWithSetMods
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
    }
}
