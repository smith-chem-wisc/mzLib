﻿using Easy.Common.Extensions;
using Omics.Fragmentation;
using System.Globalization;
using System.Text.RegularExpressions;
using Chemistry;
using Omics.Fragmentation.Peptide;
using Omics.SpectrumMatch;

namespace Readers
{
    public abstract class SpectrumMatchFromTsv
    {
        protected static readonly Regex PositionParser = new Regex(@"(\d+)\s+to\s+(\d+)");
        protected static readonly Regex VariantParser = new Regex(@"[a-zA-Z]+(\d+)([a-zA-Z]+)");
        protected static readonly Regex IonParser = new Regex(@"([a-zA-Z]+)(\d+)");

        public string FullSequence { get; protected set; }
        public int Ms2ScanNumber { get; protected set; }
        public string FileNameWithoutExtension { get; protected set; }
        public int PrecursorScanNum { get; protected set; }
        public int PrecursorCharge { get; protected set; }
        public double PrecursorMz { get; protected set; }
        public double PrecursorMass { get; protected set; }
        public double? RetentionTime { get; protected set; }
        public double Score { get; protected set; }
        public int SpectrumMatchCount { get; protected set; }
        public string Accession { get; protected set; }
        public double? SpectralAngle { get; protected set; }
        public List<MatchedFragmentIon> MatchedIons { get; protected set; }
        public Dictionary<int, List<MatchedFragmentIon>> ChildScanMatchedIons { get; protected set; }
        public double QValue { get; protected set; }
        public double PEP { get; protected set; }
        public double PEP_QValue { get; protected set; }
        public double? TotalIonCurrent { get; protected set; }
        public double? DeltaScore { get; protected set; }
        public string Notch { get; protected set; }
        public string BaseSeq { get; protected set; }
        public string EssentialSeq { get; protected set; }
        public string AmbiguityLevel { get; protected set; }
        public string MissedCleavage { get; protected set; }
        public string MonoisotopicMass { get; protected set; }
        public string MassDiffDa { get; protected set; }
        public string MassDiffPpm { get; protected set; }
        public string Name { get; protected set; }
        public string GeneName { get; protected set; }
        public string OrganismName { get; protected set; }
        public string IntersectingSequenceVariations { get; protected set; }
        public string IdentifiedSequenceVariations { get; protected set; }
        public string SpliceSites { get; protected set; }
        public string Description { get; protected set; }

        // First amino acid in protein is amino acid number 1, which differs from internal code numbering with N-terminus as 1
        // This numbering is for the peptide location within the protein
        public string StartAndEndResiduesInParentSequence { get; protected set; }
        public string PreviousResidue { get; protected set; }
        public string NextResidue { get; protected set; }
        public string DecoyContamTarget { get; protected set; }
        public double? QValueNotch { get; protected set; }

        public List<MatchedFragmentIon> VariantCrossingIons { get; protected set; }

        /// <summary>
        /// Constructor used for reading from file
        /// </summary>
        /// <param name="line"></param>
        /// <param name="split">what to split on</param>
        /// <param name="parsedHeader">index of each potential column in the header</param>
        protected SpectrumMatchFromTsv(string line, char[] split, Dictionary<string, int> parsedHeader)
        {
            var spl = line.Split(split).Select(p => p.Trim('\"')).ToArray();

            //Required properties
            FileNameWithoutExtension = spl[parsedHeader[SpectrumMatchFromTsvHeader.FileName]].Trim();

            // remove file format, e.g., .raw, .mzML, .mgf, .d
            // this is more robust but slower than Path.GetFileNameWithoutExtension
            if (FileNameWithoutExtension.Contains('.'))
            {
                foreach (var knownSpectraFileExtension in SpectrumMatchFromTsvHeader.AcceptedSpectraFormats)
                {
                    FileNameWithoutExtension = Path.GetFileName(FileNameWithoutExtension.Replace(knownSpectraFileExtension, string.Empty, StringComparison.InvariantCultureIgnoreCase));
                }
            }

            Ms2ScanNumber = int.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.Ms2ScanNumber]]);

            // this will probably not be known in an .mgf data file
            if (int.TryParse(spl[parsedHeader[SpectrumMatchFromTsvHeader.PrecursorScanNum]].Trim(), out int result))
            {
                PrecursorScanNum = result;
            }
            else
            {
                PrecursorScanNum = 0;
            }

            PrecursorCharge = (int)double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.PrecursorCharge]].Trim(), CultureInfo.InvariantCulture);
            PrecursorMz = double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.PrecursorMz]].Trim(), CultureInfo.InvariantCulture);
            PrecursorMass = double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.PrecursorMass]].Trim(), CultureInfo.InvariantCulture);
            BaseSeq = RemoveParentheses(spl[parsedHeader[SpectrumMatchFromTsvHeader.BaseSequence]].Trim());
            FullSequence = spl[parsedHeader[SpectrumMatchFromTsvHeader.FullSequence]];
            MonoisotopicMass = spl[parsedHeader[SpectrumMatchFromTsvHeader.MonoisotopicMass]].Trim();
            Score = double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.Score]].Trim(), CultureInfo.InvariantCulture);
            DecoyContamTarget = spl[parsedHeader[SpectrumMatchFromTsvHeader.DecoyContaminantTarget]].Trim();
            QValue = double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.QValue]].Trim(), CultureInfo.InvariantCulture);

            //we are reading in all primary and child ions here only to delete the child scans later. This should be done better.
            MatchedIons = (spl[parsedHeader[SpectrumMatchFromTsvHeader.MatchedIonMzRatios]].StartsWith("{")) ?
                ReadChildScanMatchedIons(spl[parsedHeader[SpectrumMatchFromTsvHeader.MatchedIonMzRatios]].Trim(), spl[parsedHeader[SpectrumMatchFromTsvHeader.MatchedIonIntensities]].Trim(), BaseSeq).First().Value :
                ReadFragmentIonsFromString(spl[parsedHeader[SpectrumMatchFromTsvHeader.MatchedIonMzRatios]].Trim(), spl[parsedHeader[SpectrumMatchFromTsvHeader.MatchedIonIntensities]].Trim(), BaseSeq, spl[parsedHeader[SpectrumMatchFromTsvHeader.MatchedIonMassDiffDa]].Trim());

            #pragma warning disable CS8601 // Possible null reference assignment.
            AmbiguityLevel = (parsedHeader[SpectrumMatchFromTsvHeader.AmbiguityLevel] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.AmbiguityLevel]].Trim();
            TotalIonCurrent = (parsedHeader[SpectrumMatchFromTsvHeader.TotalIonCurrent] < 0) ? null : (double?)double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.TotalIonCurrent]].Trim(), CultureInfo.InvariantCulture);
            DeltaScore = (parsedHeader[SpectrumMatchFromTsvHeader.DeltaScore] < 0) ? null : (double?)double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.DeltaScore]].Trim(), CultureInfo.InvariantCulture);
            Notch = (parsedHeader[SpectrumMatchFromTsvHeader.Notch] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.Notch]].Trim();
            EssentialSeq = (parsedHeader[SpectrumMatchFromTsvHeader.EssentialSequence] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.EssentialSequence]].Trim();
            MissedCleavage = (parsedHeader[SpectrumMatchFromTsvHeader.MissedCleavages] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.MissedCleavages]].Trim();
            MassDiffDa = (parsedHeader[SpectrumMatchFromTsvHeader.MassDiffDa] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.MassDiffDa]].Trim();
            MassDiffPpm = (parsedHeader[SpectrumMatchFromTsvHeader.MassDiffPpm] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.MassDiffPpm]].Trim();
            Accession = (parsedHeader[SpectrumMatchFromTsvHeader.Accession] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.Accession]].Trim();
            Name = (parsedHeader[SpectrumMatchFromTsvHeader.Name] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.Name]].Trim();
            GeneName = (parsedHeader[SpectrumMatchFromTsvHeader.GeneName] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.GeneName]].Trim();
            OrganismName = (parsedHeader[SpectrumMatchFromTsvHeader.OrganismName] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.OrganismName]].Trim();
            IntersectingSequenceVariations = (parsedHeader[SpectrumMatchFromTsvHeader.IntersectingSequenceVariations] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.IntersectingSequenceVariations]].Trim();
            IdentifiedSequenceVariations = (parsedHeader[SpectrumMatchFromTsvHeader.IdentifiedSequenceVariations] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.IdentifiedSequenceVariations]].Trim();
            SpliceSites = (parsedHeader[SpectrumMatchFromTsvHeader.SpliceSites] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.SpliceSites]].Trim();
            Description = (parsedHeader[SpectrumMatchFromTsvHeader.Description] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.Description]].Trim();
            StartAndEndResiduesInParentSequence = (parsedHeader[SpectrumMatchFromTsvHeader.StartAndEndResiduesInFullSequence] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.StartAndEndResiduesInFullSequence]].Trim();
            PreviousResidue = (parsedHeader[SpectrumMatchFromTsvHeader.PreviousResidue] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.PreviousResidue]].Trim();
            NextResidue = (parsedHeader[SpectrumMatchFromTsvHeader.NextResidue] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.NextResidue]].Trim();
            QValueNotch = (parsedHeader[SpectrumMatchFromTsvHeader.QValueNotch] < 0) ? null : (double?)double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.QValueNotch]].Trim(), CultureInfo.InvariantCulture);
            RetentionTime = (parsedHeader[SpectrumMatchFromTsvHeader.Ms2ScanRetentionTime] < 0) ? null : (double?)double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.Ms2ScanRetentionTime]].Trim(), CultureInfo.InvariantCulture);
            PEP = double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.PEP]].Trim(), CultureInfo.InvariantCulture);
            PEP_QValue = double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.PEP_QValue]].Trim(), CultureInfo.InvariantCulture);
            VariantCrossingIons = FindVariantCrossingIons();
            SpectralAngle = (parsedHeader[SpectrumMatchFromTsvHeader.SpectralAngle] < 0)
                ? null
                : (double?)double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.SpectralAngle]].Trim(),
                    CultureInfo.InvariantCulture);
            #pragma warning restore CS8601 // Possible null reference assignment.
        }

        /// <summary>
        /// Constructor used to disambiguate PsmFromTsv to a single psm object
        /// </summary>
        /// <param name="psm">psm to disambiguate</param>
        /// <param name="fullSequence">sequence of ambiguous psm to use</param>
        /// <param name="index">index of ambiguous match</param>
        protected SpectrumMatchFromTsv(SpectrumMatchFromTsv psm, string fullSequence, int index = 0, string baseSequence = "")
        {
            // psm is not ambiguous
            if (!psm.FullSequence.Contains("|"))
            {
                FullSequence = fullSequence;
                EssentialSeq = psm.EssentialSeq;
                BaseSeq = baseSequence == "" ? psm.BaseSeq : baseSequence;
                StartAndEndResiduesInParentSequence = psm.StartAndEndResiduesInParentSequence;
                Accession = psm.Accession;
                Name = psm.Name;
                GeneName = psm.GeneName;
                MonoisotopicMass = psm.MonoisotopicMass;
                MassDiffDa = psm.MassDiffDa;
                MassDiffPpm = psm.MassDiffPpm;
            }
            // potentially ambiguous fields
            else
            {
                FullSequence = fullSequence;
                EssentialSeq = psm.EssentialSeq.Split("|")[index];
                BaseSeq = baseSequence == "" ? psm.BaseSeq.Split("|")[index] : baseSequence;
                StartAndEndResiduesInParentSequence = psm.StartAndEndResiduesInParentSequence.Split("|")[index];
                Accession = psm.Accession.Split("|")[index];
                Name = psm.Name.Split("|")[index];
                GeneName = psm.GeneName.Split("|")[index];

                if (psm.MonoisotopicMass.Split("|").Count() == 1)
                {
                    MonoisotopicMass = psm.MonoisotopicMass.Split("|")[0];
                    MassDiffDa = psm.MassDiffDa.Split("|")[0];
                    MassDiffPpm = psm.MassDiffPpm.Split("|")[0];
                }
                else
                {
                    MonoisotopicMass = psm.MonoisotopicMass.Split("|")[index];
                    MassDiffDa = psm.MassDiffDa.Split("|")[index];
                    MassDiffPpm = psm.MassDiffPpm.Split("|")[index];
                }
            }

            // non ambiguous fields
            Ms2ScanNumber = psm.Ms2ScanNumber;
            FileNameWithoutExtension = psm.FileNameWithoutExtension;
            PrecursorScanNum = psm.PrecursorScanNum;
            PrecursorCharge = psm.PrecursorCharge;
            Score = psm.Score;
            MatchedIons = psm.MatchedIons.ToList();
            ChildScanMatchedIons = psm.ChildScanMatchedIons;
            QValue = psm.QValue;
            PEP = psm.PEP;
            PEP_QValue = psm.PEP_QValue;
            TotalIonCurrent = psm.TotalIonCurrent;
            DeltaScore = psm.DeltaScore;
            AmbiguityLevel = psm.AmbiguityLevel;
            MissedCleavage = psm.MissedCleavage;
            OrganismName = psm.OrganismName;
            IntersectingSequenceVariations = psm.IntersectingSequenceVariations;
            SpliceSites = psm.SpliceSites;
            Description = psm.Description;
            PreviousResidue = psm.PreviousResidue;
            NextResidue = psm.NextResidue;
            DecoyContamTarget = psm.DecoyContamTarget;
            QValueNotch = psm.QValueNotch;
            RetentionTime = psm.RetentionTime;
        }

        /// <summary>
        /// All parsing should take place within the derived class constructurs
        /// </summary>

        //Used to remove Silac labels for proper annotation
        public static string RemoveParentheses(string baseSequence)
        {
            if (baseSequence.Contains("("))
            {
                string updatedBaseSequence = "";
                bool withinParentheses = false;
                foreach (char c in baseSequence)
                {
                    if (c == ')') //leaving the parentheses
                    {
                        withinParentheses = false;
                    }
                    else if (c == '(') //entering the parentheses
                    {
                        withinParentheses = true;
                    }
                    else if (!withinParentheses) //if outside the parentheses, preserve this amino acid
                    {
                        updatedBaseSequence += c;
                    }
                    //else do nothing
                }
                return updatedBaseSequence;
            }
            return baseSequence;
        }

        /// <summary>
        /// Parses the full sequence to identify mods
        /// </summary>
        /// <param name="fullSequence"> Full sequence of the peptide in question</param>
        /// <returns> Dictionary with the key being the amino acid position of the mod and the value being the string representing the mod</returns>
        public static Dictionary<int, List<string>> ParseModifications(string fullSeq)
        {
            // use a regex to get all modifications
            string pattern = @"\[(.+?)\]";
            Regex regex = new(pattern);

            // remove each match after adding to the dict. Otherwise, getting positions
            // of the modifications will be rather difficult.
            //int patternMatches = regex.Matches(fullSeq).Count;
            Dictionary<int, List<string>> modDict = new();

            RemoveSpecialCharacters(ref fullSeq);
            MatchCollection matches = regex.Matches(fullSeq);
            int currentPosition = 0;
            foreach (Match match in matches)
            {
                GroupCollection group = match.Groups;
                string val = group[1].Value;
                int startIndex = group[0].Index;
                int captureLength = group[0].Length;
                int position = group["(.+?)"].Index;

                List<string> modList = new List<string>();
                modList.Add(val);
                // check to see if key already exist
                // if there is a missed cleavage, then there will be a label on K and a Label on X modification.
                // And, it'll be like [label]|[label] which complicates the positional stuff a little bit.
                // if the already key exists, update the current position with the capture length + 1.
                // otherwise, add the modification to the dict.

                // int to add is startIndex - current position
                int positionToAddToDict = startIndex - currentPosition;
                if (modDict.ContainsKey(positionToAddToDict))
                {
                    modDict[positionToAddToDict].Add(val);
                }
                else
                {
                    modDict.Add(positionToAddToDict, modList);
                }
                currentPosition += startIndex + captureLength;
            }
            return modDict;
        }

        /// <summary>
        /// Fixes an issue where the | appears and throws off the numbering if there are multiple mods on a single amino acid.
        /// </summary>
        /// <param name="fullSeq"></param>
        /// <param name="replacement"></param>
        /// <param name="specialCharacter"></param>
        /// <returns></returns>
        public static void RemoveSpecialCharacters(ref string fullSeq, string replacement = @"", string specialCharacter = @"\|")
        {
            // next regex is used in the event that multiple modifications are on a missed cleavage Lysine (K)
            Regex regexSpecialChar = new(specialCharacter);
            fullSeq = regexSpecialChar.Replace(fullSeq, replacement);
        }


        protected static List<MatchedFragmentIon> ReadFragmentIonsFromString(string matchedMzString, string matchedIntensityString, string peptideBaseSequence, string matchedMassErrorDaString = null, bool isProtein = true)
        {
            List<MatchedFragmentIon> matchedIons = new List<MatchedFragmentIon>();

            if (matchedMzString.Length > 2) //check if there's an ion
            {
                List<string> peakMzs = CleanMatchedIonString(matchedMzString);
                List<string> peakIntensities = CleanMatchedIonString(matchedIntensityString);
                List<string> peakMassErrorDa = null;

                if (matchedMassErrorDaString.IsNotNullOrEmpty())
                {
                    peakMassErrorDa = CleanMatchedIonString(matchedMassErrorDaString);
                }

                for (int index = 0; index < peakMzs.Count; index++)
                {
                    string peak = peakMzs[index];
                    string[] split = peak.Split(new char[] { '+', ':' }); //TODO: needs update for negative charges that doesn't break internal fragment ions or neutral losses

                    // if there is a mismatch between the number of peaks and number of intensities from the psmtsv, the intensity will be set to 1
                    double intensity = peakMzs.Count == peakIntensities.Count ? //TODO: needs update for negative charges that doesn't break internal fragment ions or neutral losses
                        double.Parse(peakIntensities[index].Split(new char[] { '+', ':', ']' })[2], CultureInfo.InvariantCulture) :
                        1.0;

                    int fragmentNumber = 0;
                    int secondaryFragmentNumber = 0;
                    ProductType productType;
                    ProductType? secondaryProductType = null;
                    FragmentationTerminus terminus = FragmentationTerminus.None; //default for internal fragments
                    int aminoAcidPosition;
                    double neutralLoss = 0;

                    //get theoretical fragment
                    string ionTypeAndNumber = split[0];

                    //if an internal fragment
                    if (ionTypeAndNumber.Contains("["))
                    {
                        // if there is no mismatch between intensity and peak counts from the psmtsv
                        if (!intensity.Equals(1.0))
                        {
                            intensity = double.Parse(peakIntensities[index].Split(new char[] { '+', ':', ']' })[3],
                                CultureInfo.InvariantCulture);
                        }
                        string[] internalSplit = split[0].Split('[');
                        string[] productSplit = internalSplit[0].Split("I");
                        string[] positionSplit = internalSplit[1].Replace("]", "").Split('-');
                        productType = (ProductType)Enum.Parse(typeof(ProductType), productSplit[0]);
                        secondaryProductType = (ProductType)Enum.Parse(typeof(ProductType), productSplit[1]);
                        fragmentNumber = int.Parse(positionSplit[0]);
                        secondaryFragmentNumber = int.Parse(positionSplit[1]);
                        aminoAcidPosition = secondaryFragmentNumber - fragmentNumber;
                    }
                    else //terminal fragment
                    {
                        Match result = IonParser.Match(ionTypeAndNumber);
                        productType = (ProductType)Enum.Parse(typeof(ProductType), result.Groups[1].Value);
                        fragmentNumber = int.Parse(result.Groups[2].Value);
                        // check for neutral loss  
                        if (ionTypeAndNumber.Contains("("))
                        {
                            string temp = ionTypeAndNumber.Replace("(", "");
                            temp = temp.Replace(")", "");
                            var split2 = temp.Split('-');
                            neutralLoss = double.Parse(split2[1], CultureInfo.InvariantCulture);
                        }

                        //get terminus
                        if (isProtein)
                            TerminusSpecificProductTypes.ProductTypeToFragmentationTerminus.TryGetValue(productType,
                                out terminus);
                        else
                            Omics.Fragmentation.Oligo.TerminusSpecificProductTypes.ProductTypeToFragmentationTerminus.TryGetValue(productType,
                                out terminus);


                        //get amino acid position
                        aminoAcidPosition = terminus is FragmentationTerminus.C or FragmentationTerminus.ThreePrime ?
                            peptideBaseSequence.Split('|')[0].Length - fragmentNumber :
                            fragmentNumber;
                    }

                    //get mass error in Daltons
                    double errorDa = 0;
                    if (matchedMassErrorDaString.IsNotNullOrEmpty() && peakMassErrorDa[index].IsNotNullOrEmpty())
                    {
                        string peakError = peakMassErrorDa[index];
                        string[] errorSplit = peakError.Split(new char[] { '+', ':', ']' });
                        errorDa = double.Parse(errorSplit[2], CultureInfo.InvariantCulture);
                    }

                    //get charge and mz
                    int z = int.Parse(split[1]);
                    double mz = double.Parse(split[2], CultureInfo.InvariantCulture);
                    double neutralExperimentalMass = mz.ToMass(z); //read in m/z converted to mass
                    double neutralTheoreticalMass = neutralExperimentalMass - errorDa; //theoretical mass is measured mass - measured error

                    //The product created here is the theoretical product, with the mass back-calculated from the measured mass and measured error
                    Product theoreticalProduct = new ProductWithCache(productType,
                      terminus,
                      neutralTheoreticalMass,
                      fragmentNumber,
                      aminoAcidPosition,
                      neutralLoss,
                      secondaryProductType,
                      secondaryFragmentNumber);

                    matchedIons.Add(new MatchedFragmentIonWithCache(theoreticalProduct, mz, intensity, z));
                }
            }
            return matchedIons;
        }
        /// <summary>
        /// Removes enclosing brackets and
        /// replaces delimimiters between ion series with comma
        /// then splits on comma
        /// </summary>
        /// <param name="input"> String containing ion series from .psmtsv </param>
        /// <returns> List of strings, with each entry containing one ion and associated property </returns>
        protected static List<string> CleanMatchedIonString(string input)
        {
            List<string> ionProperty = input.Substring(1, input.Length - 2)
                    .Replace("];[", ", ")
                    .Split(", ")
                    .ToList();
            ionProperty.RemoveAll(p => p.Contains("\"") || p.Equals(""));
            return ionProperty;
        }

        protected static Dictionary<int, List<MatchedFragmentIon>> ReadChildScanMatchedIons(string childScanMatchedMzString, string childScanMatchedIntensitiesString, string peptideBaseSequence)
        {
            var childScanMatchedIons = new Dictionary<int, List<MatchedFragmentIon>>();

            string[] matchedMzString = childScanMatchedMzString.Split(new char[] { '}' }).Where(p => !string.IsNullOrWhiteSpace(p)).ToArray();
            string[] matchedIntensityString = childScanMatchedIntensitiesString.Split(new char[] { '}' }).Where(p => !string.IsNullOrWhiteSpace(p)).ToArray();

            for (int i = 0; i < matchedMzString.Length; i++)
            {
                string[] mzsplit = matchedMzString[i].Split(new char[] { '@' });
                string[] intSplit = matchedIntensityString[i].Split(new char[] { '@' });

                int scanNumber = int.Parse(mzsplit[0].Trim(new char[] { '{' }));
                string matchedMzStrings = mzsplit[1];
                string matchedIntensityStrings = intSplit[1];

                var childMatchedIons = ReadFragmentIonsFromString(matchedMzStrings, matchedIntensityStrings, peptideBaseSequence);
                childScanMatchedIons.Add(scanNumber, childMatchedIons);
            }

            return childScanMatchedIons;
        }

        // finds the ions that contain variant residues using the position in IdentifiedSequenceVariations. When the variation spans 
        // multiple residues, if any part is contained in an ion, the ion is marked as variant crossing.
        protected List<MatchedFragmentIon> FindVariantCrossingIons()
        {
            List<MatchedFragmentIon> variantCrossingIons = new List<MatchedFragmentIon>();

            if (StartAndEndResiduesInParentSequence != null && IdentifiedSequenceVariations != null)
            {
                Match positionMatch = PositionParser.Match(StartAndEndResiduesInParentSequence);
                Match variantMatch = VariantParser.Match(IdentifiedSequenceVariations);
                if (positionMatch.Success && variantMatch.Success)
                {
                    List<ProductType> abcProductTypes = new List<ProductType>() { ProductType.a, ProductType.aDegree, ProductType.aStar,
                                                                    ProductType.b, ProductType.bWaterLoss, ProductType.bAmmoniaLoss, ProductType.c };
                    List<ProductType> xyzProductTypes = new List<ProductType>() { ProductType.x, ProductType.y, ProductType.yAmmoniaLoss,
                                                                    ProductType.yWaterLoss, ProductType.zDot, ProductType.zPlusOne};
                    int peptideStart = int.Parse(positionMatch.Groups[1].Value);
                    int peptideEnd = int.Parse(positionMatch.Groups[2].Value);
                    int variantResidueStart = int.Parse(variantMatch.Groups[1].Value);
                    int variantResidueEnd = variantResidueStart + variantMatch.Groups[2].Value.Length - 1;

                    foreach (MatchedFragmentIon ion in MatchedIons)
                    {
                        Match ionMatch = IonParser.Match(ion.Annotation);
                        if (ionMatch.Success &&
                            (variantResidueEnd >= peptideStart && variantResidueStart <= peptideEnd) &&     // variant is within peptide
                            ((abcProductTypes.Contains(ion.NeutralTheoreticalProduct.ProductType) &&        // type a, b, or c
                              peptideStart + int.Parse(ionMatch.Groups[2].Value) > variantResidueStart) ||  // crosses variant
                             (xyzProductTypes.Contains(ion.NeutralTheoreticalProduct.ProductType) &&        // type x, y, or z
                              peptideEnd - int.Parse(ionMatch.Groups[2].Value) < variantResidueEnd)))       // crosses variant
                        {
                            variantCrossingIons.Add(ion);
                        }
                    }
                }
            }
            return variantCrossingIons;
        }
        public static List<Tuple<int, string, double>> ReadLocalizedGlycan(string localizedGlycan)
        {
            List<Tuple<int, string, double>> tuples = new List<Tuple<int, string, double>>();
            if (localizedGlycan == null)
            {
                return tuples;
            }
            var lgs = localizedGlycan.Split(new string[] { "[", "]" }, StringSplitOptions.RemoveEmptyEntries);
            foreach (var lg in lgs)
            {
                var g = lg.Split(',', StringSplitOptions.RemoveEmptyEntries);

                Tuple<int, string, double> tuple = new Tuple<int, string, double>(int.Parse(g[0], CultureInfo.InvariantCulture), g[1], double.Parse(g[2], CultureInfo.InvariantCulture));
                tuples.Add(tuple);
            }

            return tuples;
        }

        public override string ToString()
        {
            return FullSequence;
        }

        public virtual LibrarySpectrum ToLibrarySpectrum()
        {
            bool isDecoy = this.DecoyContamTarget == "D";

            List<MatchedFragmentIon> fragments = new List<MatchedFragmentIon>();

            double matchedIonIntensitySum = Math.Max(1.0, this.MatchedIons.Select(i => i.Intensity).Sum());

            foreach (MatchedFragmentIon ion in this.MatchedIons)
            {
                Product product = new Product(ion.NeutralTheoreticalProduct.ProductType, ion.NeutralTheoreticalProduct.Terminus, ion.NeutralTheoreticalProduct.NeutralMass, ion.NeutralTheoreticalProduct.FragmentNumber, ion.NeutralTheoreticalProduct.AminoAcidPosition, ion.NeutralTheoreticalProduct.NeutralLoss);
                fragments.Add(new MatchedFragmentIon(product, ion.Mz, ion.Intensity / matchedIonIntensitySum, ion.Charge));
            }

            return (new(this.FullSequence, this.PrecursorMz, this.PrecursorCharge, fragments, this.RetentionTime, isDecoy));
        }
    }
}
