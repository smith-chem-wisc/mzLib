using Easy.Common.Extensions;
using Omics.Fragmentation;
using System.Globalization;
using System.Text.RegularExpressions;
using Chemistry;
using Omics.Fragmentation.Peptide;
using Omics.SpectrumMatch;
using MzLibUtil;
using System.Numerics;

namespace Readers
{
    public abstract class SpectrumMatchFromTsv : IQuantifiableRecord
    {
        protected static readonly Regex PositionParser = new Regex(@"(\d+)\s+to\s+(\d+)");
        protected static readonly Regex VariantParser = new Regex(@"[a-zA-Z]+(\d+)([a-zA-Z]+)");
        protected static readonly Regex IonParser = new Regex(@"([a-zA-Z]+)(\d+)");

        public string FullSequence { get; protected set; }
        public int Ms2ScanNumber { get; protected set; }
        public string FileNameWithoutExtension { get; protected set; }
        public int PrecursorScanNum { get; protected set; }
        public int PrecursorCharge { get; protected set; }
        public double? PrecursorIntensity { get; }
        public double PrecursorMz { get; protected set; }
        public double PrecursorMass { get; protected set; }
        public double RetentionTime { get; protected set; }
        public double Score { get; protected set; }
        public int SpectrumMatchCount { get; protected set; }
        public string Accession { get; protected set; }
        public double? SpectralAngle { get; protected set; }
        public List<MatchedFragmentIon> MatchedIons { get; set; }
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
        public string MonoisotopicMassString { get; protected set; }
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
        public double? OneOverK0 { get; protected set; }

        public List<MatchedFragmentIon> VariantCrossingIons { get; protected set; }

        #region IQuantifiableRecord Properties and Methods
        public string FileName => FileNameWithoutExtension;
        public int OneBasedScanNumber => Ms2ScanNumber;
        public string BaseSequence => BaseSeq;
        public int ChargeState => PrecursorCharge;
        public bool IsDecoy => DecoyContamTarget.Contains('D');
        public double MonoisotopicMass => double.TryParse(MonoisotopicMassString.Split('|')[0], CultureInfo.InvariantCulture, out double monoMass) ? monoMass : -1;
        private List<(string proteinAccessions, string geneName, string organism)>? _proteinGroupInfos;
        public List<(string proteinAccessions, string geneName, string organism)> ProteinGroupInfos
        {
            get
            {
                _proteinGroupInfos ??= ConstructProteinGroupInfo();
                return _proteinGroupInfos;
            }
        }
        protected List<(string proteinAccessions, string geneName, string organism)> ConstructProteinGroupInfo()
        {
            string[] accessions = Accession.Split('|');
            string[] genes = GeneName.Split('|');
            string[] organisms = OrganismName.Split('|');
            List<(string proteinAccessions, string geneName, string organism)> proteinGroupInfoList = new();
            for (int i = 0; i < accessions.Length; i++)
            {
                // Commonly, we'll have different proteins all map to the same gene and/or organism.
                // In these cases, the gene/organism column was resolved to only contain one non-ambiguous entry
                var gene = genes.Length > i ? genes[i] : genes[0];
                var organism = organisms.Length > i ? organisms[i] : organisms[0];
                proteinGroupInfoList.Add((accessions[i], gene, organism));
            }
            return proteinGroupInfoList;
        }
        

        #endregion

        /// <summary>
        /// Constructor used for reading from file
        /// </summary>
        /// <param name="line"></param>
        /// <param name="split">what to split on</param>
        /// <param name="parsedHeader">index of each potential column in the headerField</param>
        /// <param name="parsingParams">parsing parameters</param>
        protected SpectrumMatchFromTsv(string line, char[] split, Dictionary<string, int> parsedHeader, SpectrumMatchParsingParameters? parsingParams = null)
        {
            parsingParams ??= new SpectrumMatchParsingParameters();
            var spl = line.Split(split).Select(p => p.Trim('\"')).ToArray();

            //Required properties
            FileNameWithoutExtension = GetRequiredValue(SpectrumMatchFromTsvHeader.FileName, parsedHeader, spl);

            // remove file format, e.g., .raw, .mzML, .mgf, .d
            // this is more robust but slower than Path.GetFileNameWithoutExtension
            if (FileNameWithoutExtension.Contains('.'))
            {
                foreach (var knownSpectraFileExtension in SpectrumMatchFromTsvHeader.AcceptedSpectraFormats)
                {
                    FileNameWithoutExtension = Path.GetFileName(FileNameWithoutExtension.Replace(knownSpectraFileExtension, string.Empty, StringComparison.InvariantCultureIgnoreCase));
                }
            }

            Ms2ScanNumber = GetRequiredValue<int>(SpectrumMatchFromTsvHeader.Ms2ScanNumber, parsedHeader, spl);

            // this will probably not be known in an .mgf data file
            if (int.TryParse(spl[parsedHeader[SpectrumMatchFromTsvHeader.PrecursorScanNum]].Trim(), out int result))
            {
                PrecursorScanNum = result;
            }
            else
            {
                PrecursorScanNum = 0;
            }

            PrecursorCharge = (int)GetRequiredValue<double>(SpectrumMatchFromTsvHeader.PrecursorCharge, parsedHeader, spl);
            PrecursorIntensity = GetOptionalValue<double>(SpectrumMatchFromTsvHeader.PrecursorIntensity, parsedHeader, spl, null);
            PrecursorMz = GetRequiredValue<double>(SpectrumMatchFromTsvHeader.PrecursorMz, parsedHeader, spl);
            PrecursorMass = GetRequiredValue<double>(SpectrumMatchFromTsvHeader.PrecursorMass, parsedHeader, spl);
            BaseSeq = RemoveParentheses(GetRequiredValue(SpectrumMatchFromTsvHeader.BaseSequence, parsedHeader, spl));
            FullSequence = GetRequiredValue(SpectrumMatchFromTsvHeader.FullSequence, parsedHeader, spl);
            MonoisotopicMassString = GetRequiredValue(SpectrumMatchFromTsvHeader.MonoisotopicMass, parsedHeader, spl);
            Score = GetRequiredValue<double>(SpectrumMatchFromTsvHeader.Score, parsedHeader, spl);
            DecoyContamTarget = GetRequiredValue(SpectrumMatchFromTsvHeader.DecoyContaminantTarget, parsedHeader, spl);
            QValue = GetRequiredValue<double>(SpectrumMatchFromTsvHeader.QValue, parsedHeader, spl);

            //we are reading in all primary and child ions here only to delete the child scans later. This should be done better.
            MatchedIons = parsingParams.ParseMatchedFragmentIons 
                ? (spl[parsedHeader[SpectrumMatchFromTsvHeader.MatchedIonMzRatios]].StartsWith("{")) 
                    ? ReadChildScanMatchedIons(spl[parsedHeader[SpectrumMatchFromTsvHeader.MatchedIonMzRatios]].Trim(), spl[parsedHeader[SpectrumMatchFromTsvHeader.MatchedIonIntensities]].Trim(), BaseSeq, parsingParams).First().Value 
                    : ReadFragmentIonsFromString(spl[parsedHeader[SpectrumMatchFromTsvHeader.MatchedIonMzRatios]].Trim(), spl[parsedHeader[SpectrumMatchFromTsvHeader.MatchedIonIntensities]].Trim(), BaseSeq, parsingParams, spl[parsedHeader[SpectrumMatchFromTsvHeader.MatchedIonMassDiffDa]].Trim(), this is PsmFromTsv)
                : [];

            #pragma warning disable CS8601 // Possible null reference assignment.
            AmbiguityLevel = GetOptionalValue(SpectrumMatchFromTsvHeader.AmbiguityLevel, parsedHeader, spl);
            TotalIonCurrent = GetOptionalValue<double>(SpectrumMatchFromTsvHeader.TotalIonCurrent, parsedHeader, spl);
            DeltaScore = GetOptionalValue<double>(SpectrumMatchFromTsvHeader.DeltaScore, parsedHeader, spl);
            Notch = GetOptionalValue(SpectrumMatchFromTsvHeader.Notch, parsedHeader, spl);
            EssentialSeq = GetOptionalValue(SpectrumMatchFromTsvHeader.EssentialSequence, parsedHeader, spl);
            MissedCleavage = GetOptionalValue(SpectrumMatchFromTsvHeader.MissedCleavages, parsedHeader, spl);
            MassDiffDa = GetOptionalValue(SpectrumMatchFromTsvHeader.MassDiffDa, parsedHeader, spl);
            MassDiffPpm = GetOptionalValue(SpectrumMatchFromTsvHeader.MassDiffPpm, parsedHeader, spl);
            Accession = GetOptionalValue(SpectrumMatchFromTsvHeader.Accession, parsedHeader, spl);
            Name = GetOptionalValue(SpectrumMatchFromTsvHeader.Name, parsedHeader, spl);
            GeneName = GetOptionalValue(SpectrumMatchFromTsvHeader.GeneName, parsedHeader, spl);
            OrganismName = GetOptionalValue(SpectrumMatchFromTsvHeader.OrganismName, parsedHeader, spl);
            IntersectingSequenceVariations = GetOptionalValue(SpectrumMatchFromTsvHeader.IntersectingSequenceVariations, parsedHeader, spl);
            IdentifiedSequenceVariations = GetOptionalValue(SpectrumMatchFromTsvHeader.IdentifiedSequenceVariations, parsedHeader, spl);
            SpliceSites = GetOptionalValue(SpectrumMatchFromTsvHeader.SpliceSites, parsedHeader, spl);
            Description = GetOptionalValue(SpectrumMatchFromTsvHeader.Description, parsedHeader, spl);
            StartAndEndResiduesInParentSequence = GetOptionalValue(SpectrumMatchFromTsvHeader.StartAndEndResiduesInFullSequence, parsedHeader, spl);
            PreviousResidue = GetOptionalValue(SpectrumMatchFromTsvHeader.PreviousResidue, parsedHeader, spl);
            NextResidue = GetOptionalValue(SpectrumMatchFromTsvHeader.NextResidue, parsedHeader, spl);
            QValueNotch = GetOptionalValue<double>(SpectrumMatchFromTsvHeader.QValueNotch, parsedHeader, spl);
            RetentionTime = GetOptionalValue<double>(SpectrumMatchFromTsvHeader.Ms2ScanRetentionTime, parsedHeader, spl, -1).GetValueOrDefault(-1);
            PEP = GetOptionalValue<double>(SpectrumMatchFromTsvHeader.PEP, parsedHeader, spl, double.NaN).GetValueOrDefault(double.NaN);
            PEP_QValue = GetOptionalValue<double>(SpectrumMatchFromTsvHeader.PEP_QValue, parsedHeader, spl, double.NaN).GetValueOrDefault(double.NaN);
            OneOverK0 = GetOptionalValue<double>(SpectrumMatchFromTsvHeader.OneOverK0, parsedHeader, spl);
            VariantCrossingIons = FindVariantCrossingIons();
            SpectralAngle = GetOptionalValue<double>(SpectrumMatchFromTsvHeader.SpectralAngle, parsedHeader, spl);
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
                MonoisotopicMassString = psm.MonoisotopicMassString;
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

                if (psm.MonoisotopicMassString.Split("|").Count() == 1)
                {
                    MonoisotopicMassString = psm.MonoisotopicMassString.Split("|")[0];
                    MassDiffDa = psm.MassDiffDa.Split("|")[0];
                    MassDiffPpm = psm.MassDiffPpm.Split("|")[0];
                }
                else
                {
                    MonoisotopicMassString = psm.MonoisotopicMassString.Split("|")[index];
                    MassDiffDa = psm.MassDiffDa.Split("|")[index];
                    MassDiffPpm = psm.MassDiffPpm.Split("|")[index];
                }
            }

            // non ambiguous fields
            Ms2ScanNumber = psm.Ms2ScanNumber;
            FileNameWithoutExtension = psm.FileNameWithoutExtension;
            PrecursorScanNum = psm.PrecursorScanNum;
            PrecursorCharge = psm.PrecursorCharge;
            PrecursorIntensity = psm.PrecursorIntensity;
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
            OneOverK0 = psm.OneOverK0;
        }



        #region Parsing Methods

        protected string GetRequiredValue(string headerField, Dictionary<string, int> parsedHeader, string[] splitLine)
        {
            if (!parsedHeader.ContainsKey(headerField) || parsedHeader[headerField] < 0)
                throw new KeyNotFoundException($"Required Header '{headerField}' not found or invalid in parsedHeader.");

            return splitLine[parsedHeader[headerField]].Trim();
        }

        protected TNumber GetRequiredValue<TNumber>(string header, Dictionary<string, int> parsedHeader, string[] splitLine)
            where TNumber : INumber<TNumber> 
        {
            string value = GetRequiredValue(header, parsedHeader, splitLine);
            if (TNumber.TryParse(value, CultureInfo.InvariantCulture, out TNumber? result))
                return result;

            throw new FormatException($"Value '{value}' for header '{header}' could not be parsed as {typeof(TNumber).Name}.");
        }

        protected string? GetOptionalValue(string header, Dictionary<string, int> parsedHeader, string[] splitLine, string? defaultValue = null)
        {
            if (!parsedHeader.ContainsKey(header) || parsedHeader[header] < 0)
            {
                return defaultValue;
            }
            return splitLine[parsedHeader[header]].Trim();
        }

        protected TNumber? GetOptionalValue<TNumber>(string header, Dictionary<string, int> parsedHeader, string[] splitLine, TNumber? defaultValue = null)
            where TNumber : struct, IFormattable, IConvertible, IComparable, INumber<TNumber>
        {
            if (!parsedHeader.ContainsKey(header) || parsedHeader[header] < 0)
                return defaultValue;

            string value = splitLine[parsedHeader[header]].Trim();
            if (string.IsNullOrWhiteSpace(value))
                return defaultValue;

            if (TNumber.TryParse(value, CultureInfo.InvariantCulture, out TNumber result))
                return result;
            return defaultValue;
        }

        protected IHasChemicalFormula? GetOptionalChemicalFormula(string header, Dictionary<string, int> parsedHeader, string[] splitLine, IHasChemicalFormula? defaultValue = null)
        {
            if (!parsedHeader.ContainsKey(header) || parsedHeader[header] < 0)
                return defaultValue;

            string value = splitLine[parsedHeader[header]].Trim();
            if (string.IsNullOrWhiteSpace(value))
                return defaultValue;

            try
            {
                return ChemicalFormula.ParseFormula(value);
            }
            catch
            {
                return defaultValue;
            }
        }

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
        /// Parses the full sequence to identify mods. Local wrapper for MzLibUtil extension method.
        /// </summary>
        /// <param name="fullSequence"> Full sequence of the peptide in question</param>
        /// <returns> Dictionary with the key being the amino acid position of the mod and the value being the string representing the mod</returns>
        public static Dictionary<int, string> ParseModifications(string fullSeq)
        {
            return fullSeq.ParseModifications();
        }

        protected static List<MatchedFragmentIon> ReadFragmentIonsFromString(string matchedMzString, string matchedIntensityString, string peptideBaseSequence, SpectrumMatchParsingParameters parsingParams, string? matchedMassErrorDaString = null, bool isProtein = true)
        {
            List<MatchedFragmentIon> matchedIons = new List<MatchedFragmentIon>();

            if (matchedMzString.Length > 2) //check if there's an ion
            {
                List<string> peakMzs = CleanMatchedIonString(matchedMzString);
                List<string> peakIntensities = CleanMatchedIonString(matchedIntensityString);
                List<string>? peakMassErrorDa = matchedMassErrorDaString.IsNotNullOrEmpty()
                    ? CleanMatchedIonString(matchedMassErrorDaString)
                    : null;

                // Helper: Extract nth number (with sign) from a string
                static double ExtractNumber(string input, int n)
                {
                    var matches = Regex.Matches(input, @"-?\d+(\.\d+)?");
                    return matches.Count > n
                        ? double.Parse(matches[n].Value, CultureInfo.InvariantCulture)
                        : 1; // fallback default
                }

                for (int index = 0; index < peakMzs.Count; index++)
                {
                    string peak = peakMzs[index];

                    // Matches M, optional digits (to be stripped), optional custom loss, charge, m/z
                    // Examples matched: M15+1, M+1, M-P+1, M-P+1, M-A-P-H20-2, etc.
                    var mIonMatch = Regex.Match(peak, @"^(M)(\d*)([\w\-]*)([+-]\d+):([\d\.]+)$");
                    if (mIonMatch.Success)
                    {
                        // mIonMatch.Groups[1]: "M"
                        // mIonMatch.Groups[2]: digits after M (to be ignored)
                        // mIonMatch.Groups[3]: custom annotation (e.g., "-P", "-A-P-H20", or empty)
                        // mIonMatch.Groups[4]: charge (with sign)
                        // mIonMatch.Groups[5]: m/z

                        string customAnnotation = mIonMatch.Groups[3].Value; // e.g., "-P", "-A-P-H20", or ""
                        int charge = int.Parse(mIonMatch.Groups[4].Value, CultureInfo.InvariantCulture);
                        double mZ = double.Parse(mIonMatch.Groups[5].Value, CultureInfo.InvariantCulture);
                        double intens = ExtractNumber(peakIntensities[index], 2);

                        double neutralMass = mZ.ToMass(charge);
                        var product = new CustomMProduct(customAnnotation, neutralMass);
                        matchedIons.Add(parsingParams.FragmentIonsHavePlaceholderForEnvelope
                            ? new MatchedFragmentIonWithEnvelope(product, mZ, intens, charge)
                            : new MatchedFragmentIonWithCache(product, mZ, intens, charge));
                        continue;
                    }

                    // Regex: IonTypeAndNumber (with optional neutral loss), charge (+/-), m/z
                    // Examples:
                    //   y1+1:147.11267
                    //   (b5-97.98)+1:531.18657
                    //   aBaseLoss5-1:1234.489
                    //   (b5-97.98)-1:531.18657
                    var match = Regex.Match(peak, @"^(.*?)([+-]\d+):([\d\.]+)$");
                    if (!match.Success)
                        throw new FormatException($"Could not parse ion string: {peak}");

                    //get theoretical fragment
                    string ionTypeAndNumber = match.Groups[1].Value;
                    string chargeStr = match.Groups[2].Value;
                    string mzStr = match.Groups[3].Value;

                    int fragmentNumber, secondaryFragmentNumber = 0;
                    ProductType productType;
                    ProductType? secondaryProductType = null;
                    FragmentationTerminus terminus = FragmentationTerminus.None; //default for internal fragments
                    int aminoAcidPosition;
                    double neutralLoss = 0, intensity, errorDa = 0;


                    //if an internal fragment
                    if (ionTypeAndNumber.Contains('['))
                    {
                        intensity = ExtractNumber(peakIntensities[index], 3);

                        string[] internalSplit = ionTypeAndNumber.Split('[');
                        string[] productSplit = internalSplit[0].Split("I");
                        string[] positionSplit = internalSplit[1].Replace("]", "").Split('-');
                        productType = (ProductType)Enum.Parse(typeof(ProductType), productSplit[0]);
                        secondaryProductType = (ProductType)Enum.Parse(typeof(ProductType), productSplit[1]);
                        fragmentNumber = int.Parse(positionSplit[0]);
                        secondaryFragmentNumber = int.Parse(positionSplit[1]);
                        aminoAcidPosition = secondaryFragmentNumber - fragmentNumber;

                        //get mass error in Daltons
                        if (peakMassErrorDa != null && !string.IsNullOrEmpty(peakMassErrorDa[index]))
                            errorDa = ExtractNumber(peakMassErrorDa[index], 3);

                    }
                    else //terminal fragment
                    {
                        intensity = ExtractNumber(peakIntensities[index], 2);

                        Match result = IonParser.Match(ionTypeAndNumber);
                        productType = (ProductType)Enum.Parse(typeof(ProductType), result.Groups[1].Value);
                        fragmentNumber = int.Parse(result.Groups[2].Value);

                        // check for neutral loss  
                        if (ionTypeAndNumber.Contains('('))
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

                        //get mass error in Daltons
                        if (peakMassErrorDa != null && !string.IsNullOrEmpty(peakMassErrorDa[index]))
                            errorDa = ExtractNumber(peakMassErrorDa[index], 2);
                    }

                    //get charge and mz
                    int z = int.Parse(chargeStr, NumberStyles.Integer, CultureInfo.InvariantCulture);
                    double mz = double.Parse(mzStr, CultureInfo.InvariantCulture);
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

                    matchedIons.Add(parsingParams.FragmentIonsHavePlaceholderForEnvelope
                        ? new MatchedFragmentIonWithEnvelope(theoreticalProduct, mz, intensity, z)
                        : new MatchedFragmentIonWithCache(theoreticalProduct, mz, intensity, z));
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

        protected static Dictionary<int, List<MatchedFragmentIon>> ReadChildScanMatchedIons(string childScanMatchedMzString, string childScanMatchedIntensitiesString, string peptideBaseSequence, SpectrumMatchParsingParameters parsingParams)
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

                var childMatchedIons = ReadFragmentIonsFromString(matchedMzStrings, matchedIntensityStrings, peptideBaseSequence, parsingParams);
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

        #endregion

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
