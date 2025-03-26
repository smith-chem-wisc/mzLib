using Chemistry;
using Easy.Common.Extensions;
using MzLibUtil;
using Omics.Fragmentation;
using Omics.Fragmentation.Peptide;
using Omics.SpectrumMatch;
using Proteomics.PSM;
using System.Globalization;
using System.IO;
using System.Text;
using System.Text.RegularExpressions;

namespace Readers.SpectrumLibraries
{
    internal class MspSpectrumLibraryReader
    {
        private Dictionary<string, (string filePath, long byteOffset)> SequenceToFileAndLocation;
        private Dictionary<string, StreamReader> StreamReaders;
        private static Regex IonParserRegex = new Regex(@"^(\D{1,})(\d{1,})(?:[\^]|$)(\d{1,}|$)");

        /// <summary>
        /// Reads a msp file and returns LibrarySpectrum objects
        /// It is simply a cast of the ReadMsp method
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="warnings"></param>
        /// <returns></returns>
        public static List<LibrarySpectrum> ReadMspMsp(string filePath, out List<string> warnings) =>
            ReadMsp(filePath, out warnings).Cast<LibrarySpectrum>().ToList();



        /// <summary>
        /// Legacy method for reading PsmFromTsv files, creates a generic SpectrumMatchFromTsv object for each line
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="warnings"></param>
        /// <returns></returns>
        /// <exception cref="MzLibException"></exception>
        /// <exception cref="ArgumentOutOfRangeException"></exception>
        public static List<LibrarySpectrum> ReadMsp(string filePath, out List<string> warnings)
        {
            List<LibrarySpectrum> spectra = new List<LibrarySpectrum>();
            warnings = new List<string>();

            StreamReader reader = null;
            try
            {
                reader = new StreamReader(filePath);
            }
            catch (Exception e)
            {
                throw new MzLibException("Could not read file: " + e.Message, e);
            }

            int lineCount = 0;

            string line;
            Dictionary<string, int> parsedHeader = null;

            var fileType = filePath.ParseFileType();
            while (reader.Peek() > 0)
            {
                lineCount++;

                line = reader.ReadLine();

                // msp files have no header

                try
                {
                    switch (filePath.ParseFileType())
                    {
                        case SupportedFileType.msp:
                            spectra.Add(ReadLibrarySpectrum(reader));
                            break;

                        // TODO: Create an osmtsv case when transcriptomics is merged

                        default:
                            throw new ArgumentOutOfRangeException();
                    }
                }
                catch (Exception e)
                {
                    warnings.Add("Could not read line: " + lineCount);
                }
            }

            reader.Close();

            return spectra;
        }


        private static LibrarySpectrum ReadLibrarySpectrum(StreamReader reader, bool onlyReadHeader = false)
        {
            char[] nameSplit = new char[] { '/' };
            char[] mwSplit = new char[] { ':' };
            char[] commentSplit = new char[] { ' ', ':', '=' };
            char[] modSplit = new char[] { '=', '/' };
            char[] fragmentSplit = new char[] { '\t', '\"', ')', '/' };
            char[] neutralLossSplit = new char[] { '-' };

            bool readingPeaks = false;
            string sequence = null;
            int z = 2;
            double precursorMz = 0;
            double rt = 0;
            List<MatchedFragmentIon> matchedFragmentIons = new List<MatchedFragmentIon>();

            while (reader.Peek() > 0)
            {
                string line = reader.ReadLine();
                string[] split;

                if (line.StartsWith("Name", StringComparison.InvariantCultureIgnoreCase))
                {
                    //if (CrosslinkLibrarySpectrum.CrosslinkRegex.Match(line).Success)
                    //{
                    //    return ReadLibrarySpectrum_Crosslink(reader, line, onlyReadHeader);
                    //}

                    if (sequence != null)
                    {
                        return new LibrarySpectrum(sequence, precursorMz, z, matchedFragmentIons, rt);
                    }

                    split = line.Split(nameSplit);

                    // get sequence
                    sequence = split[0].Replace("Name:", string.Empty).Trim();

                    // get charge
                    z = int.Parse(split[1].Trim());
                }
                else if (line.StartsWith("MW", StringComparison.InvariantCultureIgnoreCase))
                {
                    split = line.Split(mwSplit);

                    // get precursor m/z
                    precursorMz = double.Parse(split[1].Trim(), CultureInfo.InvariantCulture);
                }
                else if (line.StartsWith("Comment", StringComparison.InvariantCultureIgnoreCase))
                {
                    split = line.Split(commentSplit);

                    // get precursor m/z if not defined yet
                    if (precursorMz == 0)
                    {
                        int indOfParent = Array.IndexOf(split, "Parent");
                        if (indOfParent > 0)
                        {
                            precursorMz = double.Parse(split[indOfParent + 1], CultureInfo.InvariantCulture);
                        }
                    }

                    // get RT
                    int indOfRt = Array.IndexOf(split, "iRT");
                    if (indOfRt > 0)
                    {
                        rt = double.Parse(split[indOfRt + 1], CultureInfo.InvariantCulture);
                    }
                    else
                    {
                        indOfRt = Array.IndexOf(split, "RT");

                        if (indOfRt > 0)
                        {
                            rt = double.Parse(split[indOfRt + 1], CultureInfo.InvariantCulture);
                        }
                    }

                    // get mods
                    // be careful about spaces! mod names can have spaces in them
                    StringBuilder sb = new StringBuilder();
                    int ind = line.IndexOf("Mods", StringComparison.InvariantCultureIgnoreCase);

                    if (ind > 0)
                    {
                        bool readingModName = false;
                        int bracketCount = 0;

                        for (int i = ind; i < line.Length; i++)
                        {
                            if (line[i] == ' ' && !readingModName)
                            {
                                break;
                            }
                            if (line[i] == '[')
                            {
                                bracketCount++;
                                readingModName = true;
                            }
                            else if (line[i] == ']')
                            {
                                bracketCount--;

                                if (bracketCount == 0)
                                {
                                    readingModName = false;
                                }
                            }

                            sb.Append(line[i]);
                        }

                        if (sb.ToString() != "Mods=0")
                        {
                            split = sb.ToString().Split(modSplit);

                            for (int i = split.Length - 1; i >= 2; i--)
                            {
                                string modString = split[i];

                                string[] modInfo = modString.Split(',');
                                int modPosition = int.Parse(modInfo[0]);
                                string modName = modInfo[2];
                                string modNameNoBrackets = modName;

                                if (modName.StartsWith('['))
                                {
                                    modNameNoBrackets = modName.Substring(1, modName.Length - 2);
                                }

                                //if (!GlobalVariables.AllModsKnownDictionary.ContainsKey(modNameNoBrackets))
                                //{
                                //    //if (PrositToMetaMorpheusModDictionary.TryGetValue(modName, out var metaMorpheusMod))
                                //    //{
                                //    //    modName = metaMorpheusMod;
                                //    //}
                                //}

                                // add the mod name into the sequence
                                string leftSeq = sequence.Substring(0, modPosition + 1);
                                string rightSeq = sequence.Substring(modPosition + 1);

                                sequence = leftSeq + modName + rightSeq;
                            }
                        }
                    }
                }
                else if (line.StartsWith("Num peaks", StringComparison.InvariantCultureIgnoreCase))
                {
                    if (onlyReadHeader)
                    {
                        return new LibrarySpectrum(sequence, precursorMz, z, matchedFragmentIons, rt);
                    }

                    // this assumes that the peaks are listed after the "Num peaks" line
                    readingPeaks = true;
                }
                else if (readingPeaks)
                {
                    matchedFragmentIons.Add(ReadFragmentIon(line, fragmentSplit, neutralLossSplit, sequence));
                }
            }

            return new LibrarySpectrum(sequence, precursorMz, z, matchedFragmentIons, rt);
        }
        /// <summary>
        /// Creates a matched fragment ion from a line in a spectral library. Does not work with P-Deep libraries.
        /// </summary>
        public static MatchedFragmentIon ReadFragmentIon(string fragmentIonLine, char[] fragmentSplit,
            char[] neutralLossSplit, string peptideSequence)
        {
            string[] split = fragmentIonLine.Split(fragmentSplit, StringSplitOptions.RemoveEmptyEntries);

            // read fragment m/z
            var experMz = double.Parse(split[0], CultureInfo.InvariantCulture);

            // read fragment intensity
            var experIntensity = double.Parse(split[1], CultureInfo.InvariantCulture);

            // read fragment type, number
            Match regexMatchResult = IonParserRegex.Match(split[2]);

            double neutralLoss = 0;
            if (split[2].Contains("-"))
            {
                String[] neutralLossInformation = split[2].Split(neutralLossSplit, StringSplitOptions.RemoveEmptyEntries).ToArray();
                neutralLoss = double.Parse(neutralLossInformation[1]);
            }

            string fragmentType = regexMatchResult.Groups[1].Value;
            int fragmentNumber = int.Parse(regexMatchResult.Groups[2].Value);
            int fragmentCharge = 1;

            if (regexMatchResult.Groups.Count > 3 && !string.IsNullOrWhiteSpace(regexMatchResult.Groups[3].Value))
            {
                fragmentCharge = int.Parse(regexMatchResult.Groups[3].Value);
            }

            ProductType peakProductType = (ProductType)Enum.Parse(typeof(ProductType), fragmentType, true);
            // Default product for productTypes not contained in the ProductTypeToFragmentationTerminus dictionary (e.g., "M" type ions)
            Product product = new Product(peakProductType, (FragmentationTerminus)Enum.Parse(typeof(FragmentationTerminus),
                "None", true), experMz, fragmentNumber, 0, 0);

            if (TerminusSpecificProductTypes.ProductTypeToFragmentationTerminus.TryGetValue(peakProductType,
                    out var terminus))
            {
                int peptideLength = peptideSequence.IsNotNullOrEmptyOrWhiteSpace() ? peptideSequence.Length : 25; // Arbitrary default peptide length
                product = new Product(peakProductType, terminus, experMz.ToMass(fragmentCharge), fragmentNumber,
                    residuePosition: terminus == FragmentationTerminus.N ? fragmentNumber : peptideLength - fragmentNumber,
                    neutralLoss);
            }

            return new MatchedFragmentIon(product, experMz, experIntensity, fragmentCharge);
        }
    }
}
