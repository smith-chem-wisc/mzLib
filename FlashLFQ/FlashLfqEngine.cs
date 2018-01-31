using Chemistry;
using MassSpectrometry;
using Proteomics;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace FlashLFQ
{
    public enum IdentificationFileType { MetaMorpheus, Morpheus, MaxQuant, TDPortal };

    public class FlashLFQEngine
    {
        // file info stuff
        public string identificationsFilePath { get; private set; }
        public string[] filePaths { get; private set; }
        public string[] analysisSummaryPerFile { get; private set; }
        public string outputFolder;

        // structures used in the FlashLFQ program
        private List<Identification> allIdentifications;
        public List<ChromatographicPeak>[] allFeaturesByFile { get; private set; }
        public HashSet<double> observedMzsToUseForIndex { get; private set; }
        public Dictionary<string, List<KeyValuePair<double, double>>> baseSequenceToIsotopicDistribution { get; private set; }
        private Dictionary<string, ProteinGroup> pgNameToProteinGroup;
        private string[] header;
        public Stopwatch globalStopwatch;
        public Stopwatch fileLocalStopwatch;

        // settings
        public bool silent { get; private set; }
        public bool pause { get; private set; }
        public int maxThreads { get; private set; }
        public IEnumerable<int> chargeStates { get; private set; }
        public double peakfindingPpmTolerance { get; private set; }
        public double ppmTolerance { get; private set; }
        public double rtTol { get; private set; }
        public double isotopePpmTolerance { get; private set; }
        public bool integrate { get; private set; }
        public int missedScansAllowed { get; private set; }
        public int numIsotopesRequired { get; private set; }
        public double mbrRtWindow { get; private set; }
        public double initialMbrRtWindow { get; private set; }
        public double mbrppmTolerance { get; private set; }
        public bool errorCheckAmbiguousMatches { get; private set; }
        public bool mbr { get; private set; }
        public bool idSpecificChargeState { get; private set; }
        public double qValueCutoff { get; private set; }
        public bool requireMonoisotopicMass { get; private set; }
        public IdentificationFileType identificationFileType { get; private set; }

        public FlashLFQEngine()
        {
            globalStopwatch = new Stopwatch();
            fileLocalStopwatch = new Stopwatch();
            allIdentifications = new List<Identification>();
            pgNameToProteinGroup = new Dictionary<string, ProteinGroup>();
            chargeStates = new List<int>();

            // default parameters
            ppmTolerance = 10.0;
            peakfindingPpmTolerance = 20.0;
            isotopePpmTolerance = 5.0;
            mbr = false;
            mbrRtWindow = 1.5;
            initialMbrRtWindow = 10.0;
            mbrppmTolerance = 5.0;
            integrate = false;
            missedScansAllowed = 1;
            numIsotopesRequired = 2;
            rtTol = 5.0;
            silent = false;
            pause = true;
            errorCheckAmbiguousMatches = true;
            maxThreads = -1;
            idSpecificChargeState = false;
            qValueCutoff = 0.01;
            requireMonoisotopicMass = true;
        }

        public bool ParseArgs(string[] args)
        {
            string[] validArgs = new string[] {
                "--idt [string|identification file path (TSV format)]",
                "--raw [string|MS data file path (.raw or .mzML allowed)]",
                "--rep [string|directory containing MS data files]",
                "--ppm [double|ppm tolerance]",
                "--iso [double|isotopic distribution tolerance in ppm]",
                "--sil [bool|silent mode]",
                "--pau [bool|pause at end of run]",
                "--int [bool|integrate features]",
                "--mbr [bool|match between runs]",
                "--chg [bool|use only precursor charge state]",
                "--rmm [bool|require observed monoisotopic mass peak]",
                "--nis [int|number of isotopes required to be observed]"
            };
            var newargs = string.Join("", args).Split(new[] { "--" }, StringSplitOptions.None);

            for (int i = 0; i < newargs.Length; i++)
                newargs[i] = newargs[i].Trim();

            newargs = newargs.Where(a => !a.Equals("")).ToArray();
            if (newargs.Length == 0)
            {
                Console.WriteLine("Accepted args are: " + string.Join(" ", validArgs) + "\n");
                Console.ReadKey();
                return false;
            }

            foreach (var arg in newargs)
            {
                try
                {
                    string flag = arg.Substring(0, 3);

                    switch (flag)
                    {
                        case ("idt"): identificationsFilePath = arg.Substring(3); break;
                        case ("raw"): filePaths = new string[] { arg.Substring(3) }; break;
                        case ("rep"):
                            string newArg = arg;
                            if (newArg.EndsWith("\""))
                                newArg = arg.Substring(0, arg.Length - 1);
                            filePaths = Directory.GetFiles(newArg.Substring(3)).Where(f => f.Substring(f.IndexOf('.')).ToUpper().Equals(".RAW") || f.Substring(f.IndexOf('.')).ToUpper().Equals(".MZML")).ToArray();
                            for (int i = 0; i < filePaths.Length; i++)
                                filePaths[i] = filePaths[i].Trim();
                            break;
                        case ("ppm"): ppmTolerance = double.Parse(arg.Substring(3)); break;
                        case ("iso"): isotopePpmTolerance = double.Parse(arg.Substring(3)); break;
                        case ("sil"): silent = Boolean.Parse(arg.Substring(3)); break;
                        case ("pau"): pause = Boolean.Parse(arg.Substring(3)); break;
                        case ("int"): integrate = Boolean.Parse(arg.Substring(3)); break;
                        case ("mbr"): mbr = Boolean.Parse(arg.Substring(3)); break;
                        case ("chg"): idSpecificChargeState = Boolean.Parse(arg.Substring(3)); break;
                        case ("rmm"): requireMonoisotopicMass = Boolean.Parse(arg.Substring(3)); break;
                        case ("nis"): numIsotopesRequired = int.Parse(arg.Substring(3)); break;
                        default:
                            if (!silent)
                            {
                                Console.WriteLine("Not a known argument: \"" + flag + "\"\n");
                                Console.WriteLine("Accepted args are: " + string.Join(" ", validArgs) + "\n");
                                Console.WriteLine("Press any key to exit");
                                Console.ReadKey();
                            }
                            return false;
                    }
                }
                catch (Exception)
                {
                    if (!silent)
                    {
                        Console.WriteLine("Can't parse argument \"" + arg + "\"\n");
                        Console.WriteLine("Accepted args are: " + string.Join(" ", validArgs) + "\n");
                        Console.WriteLine("Press any key to exit");
                        Console.ReadKey();
                    }
                    return false;
                }
            }

            if (filePaths != null && filePaths.Length == 0)
            {
                if (!silent)
                {
                    Console.WriteLine("Couldn't find any MS data files at the location specified\n");
                    Console.WriteLine("Press any key to exit");
                    Console.ReadKey();
                }
                return false;
            }

            if (outputFolder == null && identificationsFilePath != null)
                outputFolder = identificationsFilePath.Substring(0, identificationsFilePath.Length - (identificationsFilePath.Length - identificationsFilePath.IndexOf('.')));

            analysisSummaryPerFile = new string[filePaths.Length];
            allFeaturesByFile = new List<ChromatographicPeak>[filePaths.Length];
            return true;
        }

        public void PassFilePaths(string[] paths)
        {
            filePaths = paths.Distinct().ToArray();
            analysisSummaryPerFile = new string[filePaths.Length];
            allFeaturesByFile = new List<ChromatographicPeak>[filePaths.Length];
        }

        public bool ReadPeriodicTable(string optionalPeriodicTablePath)
        {
            string elementsLocation;

            try
            {
                if (optionalPeriodicTablePath == null)
                    elementsLocation = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, @"elements.dat");
                else
                    elementsLocation = optionalPeriodicTablePath;
                Loaders.LoadElements(elementsLocation);
            }
            catch (Exception)
            {
                if (!silent)
                {
                    Console.WriteLine("\nCan't read periodic table file\n");
                    Console.WriteLine("Press any key to exit");
                    Console.ReadKey();
                }
                return false;
            }
            return true;
        }

        public bool ReadIdentificationsFromTSV()
        {
            int fileNameCol = -1;
            int baseSequCol = -1;
            int fullSequCol = -1;
            int monoMassCol = -1;
            int msmsRetnCol = -1;
            int chargeStCol = -1;
            int protNameCol = -1;
            int decoyCol = -1;
            int qValueCol = -1;

            // read identification file
            if (!silent)
                Console.WriteLine("Opening identification file " + identificationsFilePath);
            string[] tsv;

            try
            {
                tsv = File.ReadAllLines(identificationsFilePath);
            }
            catch (FileNotFoundException)
            {
                if (!silent)
                {
                    Console.WriteLine("\nCan't find identification file\n");
                    Console.WriteLine("Press any key to exit");
                    Console.ReadKey();
                }
                return false;
            }
            catch (FileLoadException)
            {
                if (!silent)
                {
                    Console.WriteLine("\nCan't read identification file\n");
                    Console.WriteLine("Press any key to exit");
                    Console.ReadKey();
                }
                return false;
            }
            catch (Exception e)
            {
                if (!silent)
                {
                    Console.WriteLine("\nError reading identification file\n");
                    Console.WriteLine(e.Message + "\nPress any key to exit");
                    Console.ReadKey();
                }
                return false;
            }

            int lineCount = 0;
            foreach (var line in tsv)
            {
                lineCount++;
                var delimiters = new char[] { '\t' };

                if (lineCount == 1)
                {
                    header = line.Split(delimiters);

                    // MetaMorpheus MS/MS input
                    if (header.Contains("File Name")
                        && header.Contains("Base Sequence")
                        && header.Contains("Full Sequence")
                        && header.Contains("Peptide Monoisotopic Mass")
                        && header.Contains("Scan Retention Time")
                        && header.Contains("Precursor Charge")
                        && header.Contains("Protein Accession")
                        && header.Contains("Decoy/Contaminant/Target")
                        && header.Contains("QValue"))
                    {
                        fileNameCol = Array.IndexOf(header, "File Name");
                        baseSequCol = Array.IndexOf(header, "Base Sequence");
                        fullSequCol = Array.IndexOf(header, "Full Sequence");
                        monoMassCol = Array.IndexOf(header, "Peptide Monoisotopic Mass");
                        msmsRetnCol = Array.IndexOf(header, "Scan Retention Time");
                        chargeStCol = Array.IndexOf(header, "Precursor Charge");
                        protNameCol = Array.IndexOf(header, "Protein Accession");
                        decoyCol = Array.IndexOf(header, "Decoy/Contaminant/Target");
                        qValueCol = Array.IndexOf(header, "QValue");
                        identificationFileType = IdentificationFileType.MetaMorpheus;
                    }

                    // Morpheus MS/MS input
                    else if (header.Contains("Filename")
                        && header.Contains("Base Peptide Sequence")
                        && header.Contains("Peptide Sequence")
                        && header.Contains("Theoretical Mass (Da)")
                        && header.Contains("Retention Time (minutes)")
                        && header.Contains("Precursor Charge")
                        && header.Contains("Protein Description")
                        && header.Contains("Decoy?")
                        && header.Contains("Q-Value (%)"))
                    {
                        fileNameCol = Array.IndexOf(header, "Filename");
                        baseSequCol = Array.IndexOf(header, "Base Peptide Sequence");
                        fullSequCol = Array.IndexOf(header, "Peptide Sequence");
                        monoMassCol = Array.IndexOf(header, "Theoretical Mass (Da)");
                        msmsRetnCol = Array.IndexOf(header, "Retention Time (minutes)");
                        chargeStCol = Array.IndexOf(header, "Precursor Charge");
                        protNameCol = Array.IndexOf(header, "Protein Description");
                        decoyCol = Array.IndexOf(header, "Decoy?");
                        qValueCol = Array.IndexOf(header, "Q-Value (%)");
                        identificationFileType = IdentificationFileType.Morpheus;
                    }

                    // MaxQuant MS/MS input
                    else if (header.Contains("Raw file")
                        && header.Contains("Sequence")
                        && header.Contains("Modified sequence")
                        && header.Contains("Mass")
                        && header.Contains("Retention time")
                        && header.Contains("Charge")
                        && header.Contains("Proteins"))
                    {
                        fileNameCol = Array.IndexOf(header, "Raw file");
                        baseSequCol = Array.IndexOf(header, "Sequence");
                        fullSequCol = Array.IndexOf(header, "Modified sequence");
                        monoMassCol = Array.IndexOf(header, "Mass");
                        msmsRetnCol = Array.IndexOf(header, "Retention time");
                        chargeStCol = Array.IndexOf(header, "Charge");
                        protNameCol = Array.IndexOf(header, "Proteins");
                        identificationFileType = IdentificationFileType.MaxQuant;
                    }

                    // TDPortal Input
                    if (header.Contains("File Name")
                        && header.Contains("Sequence")
                        && header.Contains("Modifications")
                        && header.Contains("Monoisotopic Mass")
                        && header.Contains("RetentionTime")
                        && header.Contains("Accession")
                        && header.Contains("% Cleavages"))
                    {
                        fileNameCol = Array.IndexOf(header, "File Name");
                        baseSequCol = Array.IndexOf(header, "Sequence");
                        fullSequCol = Array.IndexOf(header, "Modifications");
                        monoMassCol = Array.IndexOf(header, "Monoisotopic Mass");
                        msmsRetnCol = Array.IndexOf(header, "RetentionTime");
                        protNameCol = Array.IndexOf(header, "Accession");
                        identificationFileType = IdentificationFileType.TDPortal;
                    }

                    // other search engines

                    // can't parse file
                    if (fileNameCol == -1)
                    {
                        if (!silent)
                        {
                            Console.WriteLine("Identification file is improperly formatted");
                            Console.ReadKey();
                        }
                        return false;
                    }
                }
                else
                {
                    try
                    {
                        var param = line.Split(delimiters);

                        string fileName = param[fileNameCol];
                        string BaseSequence = param[baseSequCol];
                        string ModSequence = param[fullSequCol];
                        if (identificationFileType == IdentificationFileType.TDPortal)
                            ModSequence = BaseSequence + ModSequence;
                        double monoisotopicMass = double.Parse(param[monoMassCol]);
                        double ms2RetentionTime = double.Parse(param[msmsRetnCol]);

                        int chargeState;
                        if (identificationFileType == IdentificationFileType.TDPortal)
                        {
                            chargeState = 1;
                        }
                        else
                            chargeState = (int)double.Parse(param[chargeStCol]);

                        var ident = new Identification(Path.GetFileNameWithoutExtension(fileName), BaseSequence, ModSequence, monoisotopicMass, ms2RetentionTime, chargeState);
                        allIdentifications.Add(ident);

                        ProteinGroup proteinGroup;
                        if (pgNameToProteinGroup.TryGetValue(param[protNameCol], out proteinGroup))
                            ident.proteinGroups.Add(proteinGroup);
                        else
                        {
                            proteinGroup = new ProteinGroup(param[protNameCol]);
                            pgNameToProteinGroup.Add(param[protNameCol], proteinGroup);
                            ident.proteinGroups.Add(proteinGroup);
                        }
                    }
                    catch (Exception)
                    {
                        if (!silent)
                        {
                            Console.WriteLine("Problem reading line " + lineCount + " of the TSV file");
                            Console.ReadKey();
                        }
                        return false;
                    }
                }
            }

            if (identificationFileType == IdentificationFileType.TDPortal)
            {
                var idsGroupedByFile = allIdentifications.GroupBy(p => p.fileName);
                var idFileNames = idsGroupedByFile.Select(p => p.Key);
                string[] fileNames = new string[filePaths.Length];
                for (int i = 0; i < filePaths.Length; i++)
                {
                    fileNames[i] = Path.GetFileNameWithoutExtension(filePaths[i]);
                }

                foreach (var fileName in idFileNames)
                {
                    int fileIndex = Array.IndexOf(fileNames, fileName);
                    if (fileIndex == -1)
                        continue;
                    //IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> file = OpenDataFile(fileIndex);
                    var identificationsForThisFile = allIdentifications.Where(p => Path.GetFileNameWithoutExtension(p.fileName) == Path.GetFileNameWithoutExtension(fileName));

                    foreach (var identification in identificationsForThisFile)
                    {
                        //int scanNum = file.GetClosestOneBasedSpectrumNumber(identification.ms2RetentionTime);
                        //var scan = file.GetOneBasedScan(scanNum) as IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>;
                        //if (scan != null)
                        {
                            //identification.chargeState = (int)(identification.monoisotopicMass / scan.SelectedIonMZ);
                            identification.chargeState = 30;
                        }
                    }

                    identificationsForThisFile.First().chargeState = 1;
                }
            }

            return true;
        }

        public bool Quantify(IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> file, string filePath)
        {
            if (filePaths == null)
                return false;

            fileLocalStopwatch.Restart();

            // construct bins
            var indexedMassSpectralPeaks = ConstructEmptyIndexForFile();

            // open raw file
            int i = Array.IndexOf(filePaths, filePath);
            if (i < 0)
                return false;
            var currentDataFile = file;

            // fill bins with peaks from the raw file
            var ms1ScanNumbers = IndexMassSpectralPeaks(indexedMassSpectralPeaks, currentDataFile);

            // quantify features using this file's IDs first
            allFeaturesByFile[i] = MainFileSearch(Path.GetFileNameWithoutExtension(filePath), indexedMassSpectralPeaks, ms1ScanNumbers, currentDataFile);

            // find unidentified features based on other files' identification results (MBR)
            if (mbr)
                MatchBetweenRuns(Path.GetFileNameWithoutExtension(filePath), indexedMassSpectralPeaks, allFeaturesByFile[i]);

            // error checking function
            // handles features with multiple identifying scans, and
            // also handles scans that are associated with more than one feature
            RunErrorChecking(allFeaturesByFile[i].Where(p => !p.isMbrFeature).ToList());
            allFeaturesByFile[i].RemoveAll(p => p.intensity == -1);

            foreach (var feature in allFeaturesByFile[i])
                foreach (var cluster in feature.isotopeClusters)
                    cluster.peakWithScan.Compress();

            if (!silent)
                Console.WriteLine("Finished " + Path.GetFileNameWithoutExtension(filePath));

            analysisSummaryPerFile[i] = "File analysis time = " + fileLocalStopwatch.Elapsed.ToString();
            return true;
        }

        public bool WriteResults(string baseFileName, bool writePeaks, bool writePeptides, bool writeProteins)
        {
            if (!silent)
                Console.WriteLine("Writing results");
            try
            {
                var allFeatures = allFeaturesByFile.SelectMany(p => p.Select(v => v));

                // write features
                List<string> featureOutput = new List<string> { ChromatographicPeak.TabSeparatedHeader };
                featureOutput = featureOutput.Concat(allFeatures.Select(v => v.ToString())).ToList();
                if (writePeaks)
                    File.WriteAllLines(outputFolder + baseFileName + "QuantifiedPeaks.tsv", featureOutput);

                // write baseseq groups
                var peptides = SumFeatures(allFeatures, true);
                List<string> baseSeqOutput = new List<string> { Peptide.TabSeparatedHeader };
                baseSeqOutput = baseSeqOutput.Concat(peptides.Select(p => p.ToString())).ToList();
                if (writePeptides)
                    File.WriteAllLines(outputFolder + baseFileName + "QuantifiedBaseSequences.tsv", baseSeqOutput);

                // write fullseq groups
                peptides = SumFeatures(allFeatures, false);
                List<string> fullSeqOutput = new List<string> { Peptide.TabSeparatedHeader };
                fullSeqOutput = fullSeqOutput.Concat(peptides.Select(p => p.ToString())).ToList();
                if (writePeptides)
                    File.WriteAllLines(outputFolder + baseFileName + "QuantifiedModifiedSequences.tsv", fullSeqOutput);

                // write protein results
                if (writeProteins)
                {
                    var proteinGroups = QuantifyProteins();
                    proteinGroups = proteinGroups.Where(v => v.intensitiesByFile != null).Distinct().OrderBy(p => p.proteinGroupName).ToList();
                    List<string> proteinOutput = new List<string> { string.Join("\t", new string[] { "test" }) };
                    proteinOutput = proteinOutput.Concat(proteinGroups.Select(v => v.ToString())).ToList();
                    File.WriteAllLines(outputFolder + baseFileName + "QuantifiedProteins.tsv", proteinOutput);
                }

                //write log
                List<string> logOutput = new List<string>()
                {
                    "Analysis Finish DateTime = " + DateTime.Now.ToString(),
                    "Total Analysis Time = " + globalStopwatch.Elapsed.ToString(),
                    "peakfindingPpmTolerance = " + peakfindingPpmTolerance,
                    "missedScansAllowed = " + missedScansAllowed,
                    "ppmTolerance = " + ppmTolerance,
                    "isotopePpmTolerance = " + isotopePpmTolerance,
                    "numIsotopesRequired = " + numIsotopesRequired,
                    "rtTol = " + rtTol,
                    "integrate = " + integrate,
                    "idSpecificChargeState = " + idSpecificChargeState,
                    "maxDegreesOfParallelism = " + maxThreads,
                    "mbr = " + mbr,
                    "mbrppmTolerance = " + mbrppmTolerance,
                    "mbrppmTolerance = " + mbrppmTolerance,
                    "mbrRtWindow = " + mbrRtWindow,
                    "errorCheckAmbiguousMatches = " + errorCheckAmbiguousMatches,
                    ""
                };

                for (int i = 0; i < filePaths.Length; i++)
                {
                    logOutput.Add("Analysis summary for: " + filePaths[i]);
                    logOutput.Add("\t" + analysisSummaryPerFile[i] + "\n");
                }

                File.WriteAllLines(outputFolder + baseFileName + "Log.txt", logOutput);
            }
            catch (Exception e)
            {
                if (!silent)
                {
                    Console.WriteLine("Unable to write results file to " + outputFolder);
                    Console.WriteLine(e.Message);
                    Console.WriteLine("Press any key to continue\n");
                    Console.ReadKey();
                }
                return false;
            }

            return true;
        }

        public void RetentionTimeCalibrationAndErrorCheckMatchedFeatures()
        {
            if (!silent)
                Console.WriteLine("Running retention time calibration");

            // get all unambiguous peaks for all files
            var allFeatures = allFeaturesByFile.SelectMany(p => p).Where(p => !p.isMbrFeature);
            var allAmbiguousFeatures = allFeatures.Where(p => p.numIdentificationsByFullSeq > 1).ToList();
            var ambiguousFeatureSeqs = new HashSet<string>(allAmbiguousFeatures.SelectMany(p => p.identifyingScans.Select(v => v.FullSequence)));

            foreach (var feature in allFeatures)
                if (ambiguousFeatureSeqs.Contains(feature.identifyingScans.First().FullSequence))
                    allAmbiguousFeatures.Add(feature);

            var unambiguousPeaksGroupedByFile = allFeatures.Except(allAmbiguousFeatures).Where(v => v.apexPeak != null).GroupBy(p => p.fileName);

            foreach (var file in unambiguousPeaksGroupedByFile)
            {
                // used later
                int fileIndex = 0;
                for (int i = 0; i < filePaths.Length; i++)
                {
                    if (Path.GetFileNameWithoutExtension(filePaths[i]).Equals(file.Key))
                    {
                        fileIndex = i;
                        break;
                    }
                }
                var allMbrFeaturesForThisFile = allFeaturesByFile[fileIndex].Where(p => p.isMbrFeature);

                // get the best (most intense) peak for each peptide in the file
                Dictionary<string, ChromatographicPeak> pepToBestFeatureForThisFile = new Dictionary<string, ChromatographicPeak>();
                foreach (var testPeak in file)
                {
                    ChromatographicPeak currentBestPeak;
                    if (pepToBestFeatureForThisFile.TryGetValue(testPeak.identifyingScans.First().FullSequence, out currentBestPeak))
                    {
                        if (currentBestPeak.intensity > testPeak.intensity)
                            pepToBestFeatureForThisFile[testPeak.identifyingScans.First().FullSequence] = testPeak;
                    }
                    else
                        pepToBestFeatureForThisFile.Add(testPeak.identifyingScans.First().FullSequence, testPeak);
                }


                foreach (var otherFile in unambiguousPeaksGroupedByFile)
                {
                    // get the other files' best peak for the same peptides (to make an RT calibration curve) 
                    if (otherFile.Key.Equals(file.Key))
                        continue;

                    var featuresInCommon = otherFile.Where(p => pepToBestFeatureForThisFile.ContainsKey(p.identifyingScans.First().FullSequence));

                    Dictionary<string, ChromatographicPeak> pepToBestFeatureForOtherFile = new Dictionary<string, ChromatographicPeak>();
                    foreach (var testPeak in featuresInCommon)
                    {
                        ChromatographicPeak currentBestPeak;
                        if (pepToBestFeatureForOtherFile.TryGetValue(testPeak.identifyingScans.First().FullSequence, out currentBestPeak))
                        {
                            if (currentBestPeak.intensity > testPeak.intensity)
                                pepToBestFeatureForOtherFile[testPeak.identifyingScans.First().FullSequence] = testPeak;
                        }
                        else
                            pepToBestFeatureForOtherFile.Add(testPeak.identifyingScans.First().FullSequence, testPeak);
                    }

                    // create a rt-to-rt correlation for the two files' peptides
                    Dictionary<string, Tuple<double, double>> rtCalPoints = new Dictionary<string, Tuple<double, double>>();

                    foreach (var kvp in pepToBestFeatureForOtherFile)
                        rtCalPoints.Add(kvp.Key, new Tuple<double, double>(pepToBestFeatureForThisFile[kvp.Key].apexPeak.peakWithScan.retentionTime, kvp.Value.apexPeak.peakWithScan.retentionTime));

                    var someDoubles = rtCalPoints.Select(p => (p.Value.Item1 - p.Value.Item2));
                    double average = someDoubles.Average();
                    double sumOfSquaresOfDifferences = someDoubles.Select(val => (val - average) * (val - average)).Sum();
                    double sd = Math.Sqrt(sumOfSquaresOfDifferences / (someDoubles.Count() - 1));


                    // remove extreme outliers
                    if (sd > 1.0)
                    {
                        var pointsToRemove = rtCalPoints.Where(p => p.Value.Item1 - p.Value.Item2 > average + sd || p.Value.Item1 - p.Value.Item2 < average - sd).ToList();
                        foreach (var point in pointsToRemove)
                            rtCalPoints.Remove(point.Key);

                        someDoubles = rtCalPoints.Select(p => (p.Value.Item1 - p.Value.Item2));
                        average = someDoubles.Average();
                        sumOfSquaresOfDifferences = someDoubles.Select(val => (val - average) * (val - average)).Sum();
                        sd = Math.Sqrt(sumOfSquaresOfDifferences / (someDoubles.Count() - 1));
                    }


                    // find rt differences between files
                    List<Tuple<double, double>> rtCalPoints2 = rtCalPoints.Values.OrderBy(p => p.Item1).ToList();

                    int minRt = (int)Math.Floor(rtCalPoints2.First().Item1);
                    int maxRt = (int)Math.Ceiling(rtCalPoints2.Last().Item1);
                    var roughRts = Enumerable.Range(minRt, (maxRt - minRt) + 1);
                    var rtCalibrationDictionary = new Dictionary<int, List<double>>();
                    foreach (var rt in rtCalPoints2)
                    {
                        List<double> points = null;
                        int intRt = (int)Math.Round(rt.Item1);
                        if (rtCalibrationDictionary.TryGetValue(intRt, out points))
                            points.Add(rt.Item1 - rt.Item2);
                        else
                            rtCalibrationDictionary.Add(intRt, new List<double> { rt.Item1 - rt.Item2 });
                    }

                    List<string> output = new List<string>();
                    foreach (var point in rtCalPoints)
                    {
                        output.Add("" + point.Key + "\t" + point.Value.Item1 + "\t" + point.Value.Item2 + "\t" + (point.Value.Item1 - point.Value.Item2));
                    }

                    //File.WriteAllLines(outputFolder + file.Key + otherFile.Key + "RTCal.tsv", output);

                    output = new List<string>();
                    // create spline
                    // index is minute of source, double is rt calibration factor (in minutes) for destination file
                    double[] rtCalRunningSpline = new double[maxRt + 1];
                    double[] stdevRunningSpline = new double[maxRt + 1];
                    for (int i = 0; i < rtCalRunningSpline.Length; i++)
                    {
                        rtCalRunningSpline[i] = double.NaN;
                        stdevRunningSpline[i] = double.NaN;
                    }

                    // calculate stdev for each element in spline
                    for (int i = 1; i <= maxRt; i++)
                    {
                        List<double> list;
                        if (rtCalibrationDictionary.TryGetValue(i, out list))
                        {
                            if (list.Count > 3)
                            {
                                list.Sort();
                                rtCalRunningSpline[i] = list[list.Count / 2]; // median rt difference for this timepoint

                                average = list.Average();
                                sumOfSquaresOfDifferences = list.Select(val => (val - average) * (val - average)).Sum();
                                sd = Math.Sqrt(sumOfSquaresOfDifferences / (list.Count - 1));

                                if (3 * sd > (mbrRtWindow / 2.0))
                                    stdevRunningSpline[i] = mbrRtWindow / 2.0;
                                else
                                    stdevRunningSpline[i] = 3 * sd;
                            }
                        }
                    }

                    // fill gaps in spline (linear interpolation)
                    for (int i = minRt; i <= maxRt; i++)
                    {
                        if (double.IsNaN(rtCalRunningSpline[i]))
                        {
                            KeyValuePair<int, double> prevCalPoint = new KeyValuePair<int, double>(0, double.NaN);
                            KeyValuePair<int, double> nextCalPoint = new KeyValuePair<int, double>(0, double.NaN);

                            for (int j = i; j >= minRt; j--)
                            {
                                if (j < i - 3)
                                    break;
                                if (!double.IsNaN(rtCalRunningSpline[j]))
                                {
                                    prevCalPoint = new KeyValuePair<int, double>(j, rtCalRunningSpline[j]);
                                    break;
                                }
                            }
                            for (int j = i; j <= maxRt; j++)
                            {
                                if (j > i + 3)
                                    break;
                                if (!double.IsNaN(rtCalRunningSpline[j]))
                                {
                                    nextCalPoint = new KeyValuePair<int, double>(j, rtCalRunningSpline[j]);
                                    break;
                                }
                            }

                            if (!double.IsNaN(prevCalPoint.Value) && !double.IsNaN(nextCalPoint.Value))
                            {
                                var slope = (prevCalPoint.Value - nextCalPoint.Value) / (prevCalPoint.Key - nextCalPoint.Key);
                                var yint = prevCalPoint.Value - slope * prevCalPoint.Key;
                                double interpolatedRtCalPoint = slope * i + yint;

                                rtCalRunningSpline[i] = interpolatedRtCalPoint;
                                stdevRunningSpline[i] = stdevRunningSpline[i] = mbrRtWindow / 2.0;
                            }
                        }
                    }

                    for (int i = minRt; i <= maxRt; i++)
                        output.Add("" + i + "\t" + rtCalRunningSpline[i] + "\t" + stdevRunningSpline[i]);

                    //File.WriteAllLines(outputFolder + file.Key + otherFile.Key + "RTCal2.tsv", output);


                    // finished rt calibration for these 2 files; now use rt cal spline to find matched features
                    var allMatchedFeaturesToLookForNow = allMbrFeaturesForThisFile.Where(p => p.identifyingScans.First().fileName.Equals(otherFile.Key)).ToList();

                    // filter peak candidates with rt cal to get apex peak
                    foreach (var mbrFeature in allMatchedFeaturesToLookForNow)
                    {
                        if (mbrFeature.isotopeClusters.Any())
                        {
                            // shift = thisFileRt - otherFileRt
                            int rtSplineLookupTime = (int)Math.Round(mbrFeature.identifyingScans.First().ms2RetentionTime);

                            if (rtSplineLookupTime < rtCalRunningSpline.Length)
                            {
                                double rtShift = rtCalRunningSpline[rtSplineLookupTime];
                                double rtToleranceHere = stdevRunningSpline[rtSplineLookupTime];
                                double theoreticalRt = mbrFeature.identifyingScans.First().ms2RetentionTime + rtShift;

                                if (!double.IsNaN(rtShift))
                                    mbrFeature.isotopeClusters = mbrFeature.isotopeClusters.Where(p => Math.Abs(p.peakWithScan.retentionTime - theoreticalRt) < rtToleranceHere).ToList();
                                else
                                    mbrFeature.isotopeClusters = new List<IsotopeCluster>();
                            }
                            else
                                mbrFeature.isotopeClusters = new List<IsotopeCluster>();
                        }
                    }

                    output = new List<string>();
                    foreach (var feature in allMatchedFeaturesToLookForNow)
                    {
                        if (feature.isotopeClusters.Any())
                        {
                            feature.CalculateIntensityForThisFeature(integrate);
                            output.Add("" + feature.apexPeak.peakWithScan.retentionTime + "\t" + feature.identifyingScans.First().ms2RetentionTime);
                        }
                        else
                            feature.intensity = -1;
                    }

                    //File.WriteAllLines(outputFolder + file.Key + otherFile.Key + "RTCal3.tsv", output);
                }
            }

            for (int i = 0; i < allFeaturesByFile.Length; i++)
            {
                // remove empty mbr features
                allFeaturesByFile[i].RemoveAll(p => p.intensity == -1);
                allFeaturesByFile[i] = allFeaturesByFile[i].Where(p => !p.isMbrFeature || (p.isMbrFeature && p.isotopeClusters.Any())).ToList();
                // error check
                RunErrorChecking(allFeaturesByFile[i]);
            }
        }

        public List<ProteinGroup> QuantifyProteins()
        {
            List<ProteinGroup> returnList = new List<ProteinGroup>();

            if (!silent)
                Console.WriteLine("Quantifying proteins");
            var fileNames = filePaths.Select(p => Path.GetFileNameWithoutExtension(p)).ToList();

            var allFeatures = allFeaturesByFile.SelectMany(p => p);
            var allAmbiguousFeatures = allFeatures.Where(p => p.numIdentificationsByBaseSeq > 1).ToList();
            var ambiguousFeatureSeqs = new HashSet<string>(allAmbiguousFeatures.SelectMany(p => p.identifyingScans.Select(v => v.BaseSequence)));

            foreach (var feature in allFeatures)
            {
                if (ambiguousFeatureSeqs.Contains(feature.identifyingScans.First().BaseSequence))
                    allAmbiguousFeatures.Add(feature);
            }

            var allUnambiguousFeatures = allFeatures.Except(allAmbiguousFeatures).ToList();

            var proteinsWithFeatures = new Dictionary<ProteinGroup, List<ChromatographicPeak>>();
            foreach(var feature in allUnambiguousFeatures)
            {
                foreach(var proteinGroup in feature.identifyingScans.First().proteinGroups)
                {
                    List<ChromatographicPeak> featuresForThisProtein;

                    if (proteinsWithFeatures.TryGetValue(proteinGroup, out featuresForThisProtein))
                        featuresForThisProtein.Add(feature);
                    else
                        proteinsWithFeatures.Add(proteinGroup, new List<ChromatographicPeak> { feature });
                }
            }

            foreach(var proteinGroupsFeatures in proteinsWithFeatures)
            {
                var peaksForThisProteinPerFile = proteinGroupsFeatures.Value.GroupBy(p => p.fileName);
                proteinGroupsFeatures.Key.intensitiesByFile = new double[fileNames.Count];

                foreach (var file in peaksForThisProteinPerFile)
                {
                    int i = fileNames.IndexOf(file.Key);
                    foreach (var feature in file)
                    {
                        int numProteinGroupsClaimingThisFeature = feature.identifyingScans.SelectMany(p => p.proteinGroups).Distinct().Count();
                        proteinGroupsFeatures.Key.intensitiesByFile[i] += (feature.intensity / numProteinGroupsClaimingThisFeature);
                    }
                }

                returnList.Add(proteinGroupsFeatures.Key);
            }

            return returnList;
        }

        public void AddIdentification(string fileName, string BaseSequence, string FullSequence, double monoisotopicMass, double ms2RetentionTime, int chargeState, List<string> proteinGroupNames)
        {
            var ident = new Identification(fileName, BaseSequence, FullSequence, monoisotopicMass, ms2RetentionTime, chargeState);

            foreach (var pgGroupName in proteinGroupNames)
            {
                ProteinGroup pg;
                if (pgNameToProteinGroup.TryGetValue(pgGroupName, out pg))
                    ident.proteinGroups.Add(pg);
                else
                {
                    pg = new ProteinGroup(pgGroupName);
                    pgNameToProteinGroup.Add(pgGroupName, pg);
                    ident.proteinGroups.Add(pg);
                }
            }

            allIdentifications.Add(ident);
        }
        
        public void ConstructIndexTemplateFromIdentifications()
        {
            observedMzsToUseForIndex = new HashSet<double>();

            var peptideGroups = allIdentifications.GroupBy(p => p.FullSequence).ToList();
            var peptideBaseSeqs = new HashSet<string>(allIdentifications.Select(p => p.BaseSequence));
            var numCarbonsToIsotopicDistribution = new Dictionary<int, IsotopicDistribution>();
            baseSequenceToIsotopicDistribution = new Dictionary<string, List<KeyValuePair<double, double>>>();

            foreach (var baseSeq in peptideBaseSeqs)
            {
                if (baseSequenceToIsotopicDistribution.ContainsKey(baseSeq))
                    continue;

                Proteomics.Peptide p = new Proteomics.Peptide(baseSeq);
                int numCarbonsInThisPeptide = p.ElementCountWithIsotopes("C");

                // get expected C13 mass shifts and abundances
                IsotopicDistribution isotopicDistribution;
                if (!numCarbonsToIsotopicDistribution.TryGetValue(numCarbonsInThisPeptide, out isotopicDistribution))
                {
                    isotopicDistribution = IsotopicDistribution.GetDistribution(ChemicalFormula.ParseFormula("C" + numCarbonsInThisPeptide), 0.00001, 0.001);
                    numCarbonsToIsotopicDistribution.Add(numCarbonsInThisPeptide, isotopicDistribution);
                }

                var masses = isotopicDistribution.Masses.ToArray();
                var abundances = isotopicDistribution.Intensities.ToArray();
                var isotopicMassesAndNormalizedAbundances = new List<KeyValuePair<double, double>>();

                var monoisotopicMass = masses.Min();
                var highestAbundance = abundances.Max();

                for (int i = 0; i < masses.Length; i++)
                {
                    // expected isotopic mass shifts for peptide of this length
                    masses[i] -= monoisotopicMass;

                    // normalized abundance of each mass shift
                    abundances[i] /= highestAbundance;

                    // look for these isotopes
                    if (i < (numIsotopesRequired - 1) || abundances[i] > 0.1)
                        isotopicMassesAndNormalizedAbundances.Add(new KeyValuePair<double, double>(masses[i], abundances[i]));
                }

                baseSequenceToIsotopicDistribution.Add(baseSeq, isotopicMassesAndNormalizedAbundances);
            }

            var minChargeState = allIdentifications.Select(p => p.chargeState).Min();
            var maxChargeState = allIdentifications.Select(p => p.chargeState).Max();
            chargeStates = Enumerable.Range(minChargeState, (maxChargeState - minChargeState) + 1);

            // build theoretical m/z bins
            foreach (var pepGroup in peptideGroups)
            {
                double lowestCommonMassShift = baseSequenceToIsotopicDistribution[pepGroup.First().BaseSequence].Select(p => p.Key).Min();
                var mostCommonIsotopeShift = baseSequenceToIsotopicDistribution[pepGroup.First().BaseSequence].Where(p => p.Value == 1).First().Key;

                var thisPeptidesLowestCommonMass = pepGroup.First().monoisotopicMass + lowestCommonMassShift;
                var thisPeptidesMostAbundantMass = pepGroup.First().monoisotopicMass + mostCommonIsotopeShift;

                foreach (var pep in pepGroup)
                {
                    pep.massToLookFor = requireMonoisotopicMass ? pepGroup.First().monoisotopicMass : thisPeptidesMostAbundantMass;
                }

                foreach (var chargeState in chargeStates)
                {
                    var t = ClassExtensions.ToMz(pepGroup.First().massToLookFor, chargeState);
                    double floorMz = Math.Floor(t * 100) / 100;
                    double ceilingMz = Math.Ceiling(t * 100) / 100;

                    if (!observedMzsToUseForIndex.Contains(floorMz))
                        observedMzsToUseForIndex.Add(floorMz);
                    if (!observedMzsToUseForIndex.Contains(ceilingMz))
                        observedMzsToUseForIndex.Add(ceilingMz);
                }
            }
        }

        public void SetParallelization(int maxDegreesOfParallelism)
        {
            this.maxThreads = maxDegreesOfParallelism;
        }

        private Dictionary<double, List<IndexedMassSpectralPeak>> ConstructEmptyIndexForFile()
        {
            return observedMzsToUseForIndex.ToDictionary(v => v, v => new List<IndexedMassSpectralPeak>());
        }

        private List<KeyValuePair<int, double>> IndexMassSpectralPeaks(Dictionary<double, List<IndexedMassSpectralPeak>> mzBins, IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> file)
        {
            var ms1ScanNumbersWithRetentionTimes = new List<KeyValuePair<int, double>>();

            var allMs1Scans = file.GetMS1Scans().ToList();

            foreach (var scan in allMs1Scans)
                ms1ScanNumbersWithRetentionTimes.Add(new KeyValuePair<int, double>(scan.OneBasedScanNumber, scan.RetentionTime));
            ms1ScanNumbersWithRetentionTimes = ms1ScanNumbersWithRetentionTimes.OrderBy(p => p.Key).ToList();

            if (!silent)
                Console.WriteLine("Assigning MS1 peaks to bins");

            //multithreaded bin-filling
            var allGoodPeaks = new List<List<KeyValuePair<double, IndexedMassSpectralPeak>>>();

            Parallel.ForEach(Partitioner.Create(0, allMs1Scans.Count),
                new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                (range, loopState) =>
                {
                    var threadLocalGoodPeaks = new List<KeyValuePair<double, IndexedMassSpectralPeak>>();

                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        int peakIndexInThisScan = 0;

                        for (int j = 0; j < allMs1Scans[i].MassSpectrum.XArray.Length; j++)
                        {
                            IndexedMassSpectralPeak element = null;
                            double floorMz = (Math.Floor(allMs1Scans[i].MassSpectrum.XArray[j] * 100) / 100);
                            double ceilingMz = (Math.Ceiling(allMs1Scans[i].MassSpectrum.XArray[j] * 100) / 100);

                            if (mzBins.ContainsKey(floorMz))
                            {
                                element = new IndexedMassSpectralPeak(new MassSpectralPeak(allMs1Scans[i].MassSpectrum.XArray[j], allMs1Scans[i].MassSpectrum.YArray[j]), allMs1Scans[i], peakIndexInThisScan);
                                threadLocalGoodPeaks.Add(new KeyValuePair<double, IndexedMassSpectralPeak>(floorMz, element));
                            }
                            if (mzBins.ContainsKey(ceilingMz))
                            {
                                if (element == null)
                                    element = new IndexedMassSpectralPeak(new MassSpectralPeak(allMs1Scans[i].MassSpectrum.XArray[j], allMs1Scans[i].MassSpectrum.YArray[j]), allMs1Scans[i], peakIndexInThisScan);
                                threadLocalGoodPeaks.Add(new KeyValuePair<double, IndexedMassSpectralPeak>(ceilingMz, element));
                            }

                            peakIndexInThisScan++;
                        }
                    }

                    lock (allGoodPeaks)
                        allGoodPeaks.Add(threadLocalGoodPeaks);
                }
            );

            Parallel.ForEach(Partitioner.Create(0, allGoodPeaks.Count),
                new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                (range) =>
                {
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        foreach (var element in allGoodPeaks[i])
                        {
                            var t = mzBins[element.Key];
                            lock (t)
                                t.Add(element.Value);
                        }
                    }
                });

            /*
            // single-threaded bin-filling
            foreach (var scan in allMs1Scans)
            {
                int peakIndexInThisScan = 0;

                foreach (var peak in scan.MassSpectrum)
                {
                    List<FlashLFQMzBinElement> mzBin;
                    FlashLFQMzBinElement element = null;
                    double floorMz = Math.Floor(peak.Mz * 100) / 100;
                    double ceilingMz = Math.Ceiling(peak.Mz * 100) / 100;

                    if (mzBins.TryGetValue(floorMz, out mzBin))
                    {
                        element = new FlashLFQMzBinElement(peak, scan, peakIndexInThisScan);
                        mzBin.Add(element);
                    }

                    if (mzBins.TryGetValue(ceilingMz, out mzBin))
                    {
                        if (element == null)
                            element = new FlashLFQMzBinElement(peak, scan, peakIndexInThisScan);
                        mzBin.Add(element);
                    }

                    peakIndexInThisScan++;
                }
            }
            */

            return ms1ScanNumbersWithRetentionTimes;
        }

        private List<ChromatographicPeak> MainFileSearch(string fileName, Dictionary<double, List<IndexedMassSpectralPeak>> mzBins, List<KeyValuePair<int, double>> ms1ScanNumbersWithRts, IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> file)
        {
            if (!silent)
                Console.WriteLine("Quantifying peptides for " + fileName);

            var ms1ScanNumbers = ms1ScanNumbersWithRts.Select(p => p.Key).OrderBy(p => p).ToList();
            var concurrentBagOfFeatures = new ConcurrentBag<ChromatographicPeak>();

            var groups = allIdentifications.GroupBy(p => p.fileName);
            var identificationsForThisFile = groups.Where(p => p.Key.Equals(fileName)).FirstOrDefault();

            if (identificationsForThisFile == null)
                return concurrentBagOfFeatures.ToList();

            var identifications = identificationsForThisFile.ToList();

            Parallel.ForEach(Partitioner.Create(0, identifications.Count),
                new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                (range, loopState) =>
                {
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        var identification = identifications[i];
                        ChromatographicPeak msmsFeature = new ChromatographicPeak();
                        msmsFeature.identifyingScans.Add(identification);
                        msmsFeature.isMbrFeature = false;
                        msmsFeature.fileName = fileName;

                        foreach (var chargeState in chargeStates)
                        {
                            if (idSpecificChargeState)
                                if (chargeState != identification.chargeState)
                                    continue;

                            double theorMzHere = ClassExtensions.ToMz(identification.massToLookFor, chargeState);
                            double mzTolHere = ((peakfindingPpmTolerance / 1e6) * identification.monoisotopicMass) / chargeState;

                            double floorMz = Math.Floor(theorMzHere * 100) / 100;
                            double ceilingMz = Math.Ceiling(theorMzHere * 100) / 100;

                            IEnumerable<IndexedMassSpectralPeak> binPeaks = new List<IndexedMassSpectralPeak>();
                            List<IndexedMassSpectralPeak> list;

                            for (double j = floorMz; j <= ceilingMz; j += 0.01)
                            {
                                if (mzBins.TryGetValue(Math.Round(j, 2), out list))
                                    binPeaks = binPeaks.Concat(list);
                            }

                            // filter by mz tolerance
                            var binPeaksHere = binPeaks.Where(p => Math.Abs(p.mainPeak.Mz - theorMzHere) < mzTolHere);
                            // remove duplicates
                            binPeaksHere = binPeaksHere.Distinct();
                            // filter by RT
                            binPeaksHere = binPeaksHere.Where(p => Math.Abs(p.retentionTime - identification.ms2RetentionTime) < rtTol);

                            if (binPeaksHere.Any())
                            {
                                // get precursor scan to start at
                                int precursorScanNum = 0;
                                foreach (var ms1Scan in ms1ScanNumbersWithRts)
                                {
                                    if (ms1Scan.Value < identification.ms2RetentionTime)
                                        precursorScanNum = ms1Scan.Key;
                                    else
                                        break;
                                }
                                if (precursorScanNum == 0)
                                    throw new Exception("Error getting precursor scan number");

                                // separate peaks by rt into left and right of the identification RT
                                var rightPeaks = binPeaksHere.Where(p => p.retentionTime >= identification.ms2RetentionTime).OrderBy(p => p.retentionTime);
                                var leftPeaks = binPeaksHere.Where(p => p.retentionTime < identification.ms2RetentionTime).OrderByDescending(p => p.retentionTime);

                                // store peaks on each side of the identification RT
                                var crawledRightPeaks = ScanCrawl(rightPeaks, missedScansAllowed, precursorScanNum, ms1ScanNumbers);
                                var crawledLeftPeaks = ScanCrawl(leftPeaks, missedScansAllowed, precursorScanNum, ms1ScanNumbers);

                                // filter again by smaller mz tolerance
                                mzTolHere = ((ppmTolerance / 1e6) * identification.monoisotopicMass) / chargeState;
                                var validPeaks = crawledRightPeaks.Concat(crawledLeftPeaks);
                                validPeaks = validPeaks.Where(p => Math.Abs(p.mainPeak.Mz - theorMzHere) < mzTolHere);

                                // filter by isotopic distribution
                                var validIsotopeClusters = FilterPeaksByIsotopicDistribution(validPeaks, identification, chargeState, true);

                                // if multiple mass spectral peaks in the same scan are valid, pick the one with the smallest mass error
                                var peaksInSameScan = validIsotopeClusters.GroupBy(p => p.peakWithScan.oneBasedScanNumber).Where(v => v.Count() > 1);
                                if (peaksInSameScan.Any())
                                {
                                    foreach (var group in peaksInSameScan)
                                    {
                                        var mzToUse = group.Select(p => Math.Abs(p.peakWithScan.mainPeak.Mz - theorMzHere)).Min();
                                        var peakToUse = group.Where(p => Math.Abs(p.peakWithScan.mainPeak.Mz - theorMzHere) == mzToUse).First();
                                        var peaksToRemove = group.Where(p => p != peakToUse);
                                        validIsotopeClusters = validIsotopeClusters.Except(peaksToRemove);
                                    }
                                }

                                foreach (var validCluster in validIsotopeClusters)
                                    msmsFeature.isotopeClusters.Add(validCluster);
                            }
                        }

                        msmsFeature.CalculateIntensityForThisFeature(integrate);
                        CutPeak(msmsFeature, integrate, ms1ScanNumbers);
                        concurrentBagOfFeatures.Add(msmsFeature);
                    }
                }
            );

            // merge results from all threads together
            return concurrentBagOfFeatures.ToList();
        }

        private void MatchBetweenRuns(string thisFileName, Dictionary<double, List<IndexedMassSpectralPeak>> mzBins, List<ChromatographicPeak> features)
        {
            if (!silent)
                Console.WriteLine("Finding possible matched peptides for " + thisFileName);

            var concurrentBagOfMatchedFeatures = new ConcurrentBag<ChromatographicPeak>();
            var identificationsFromOtherRunsToLookFor = new List<Identification>();
            var idsGroupedByFullSeq = allIdentifications.GroupBy(p => p.FullSequence);

            foreach (var fullSequenceGroup in idsGroupedByFullSeq)
            {
                // look for peptides with no ID's in this file
                var seqsByFilename = fullSequenceGroup.GroupBy(p => p.fileName);

                if (!seqsByFilename.Where(p => p.Key.Equals(thisFileName)).Any())
                    identificationsFromOtherRunsToLookFor.AddRange(fullSequenceGroup);
            }

            if (identificationsFromOtherRunsToLookFor.Any())
            {
                Parallel.ForEach(Partitioner.Create(0, identificationsFromOtherRunsToLookFor.Count),
                    new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                    (range, loopState) =>
                    {
                        for (int i = range.Item1; i < range.Item2; i++)
                        {
                            var identification = identificationsFromOtherRunsToLookFor[i];

                            ChromatographicPeak mbrFeature = new ChromatographicPeak();
                            mbrFeature.identifyingScans.Add(identification);
                            mbrFeature.isMbrFeature = true;
                            mbrFeature.fileName = thisFileName;

                            foreach (var chargeState in chargeStates)
                            {
                                double theorMzHere = ClassExtensions.ToMz(identification.massToLookFor, chargeState);
                                double mzTolHere = ((mbrppmTolerance / 1e6) * identification.monoisotopicMass) / chargeState;

                                double floorMz = Math.Floor(theorMzHere * 100) / 100;
                                double ceilingMz = Math.Ceiling(theorMzHere * 100) / 100;

                                IEnumerable<IndexedMassSpectralPeak> binPeaks = new List<IndexedMassSpectralPeak>();
                                List<IndexedMassSpectralPeak> list;
                                for (double j = floorMz; j <= ceilingMz; j += 0.01)
                                {
                                    if (mzBins.TryGetValue(Math.Round(j, 2), out list))
                                        binPeaks = binPeaks.Concat(list);
                                }

                                // filter by mz tolerance
                                var binPeaksHere = binPeaks.Where(p => Math.Abs(p.mainPeak.Mz - theorMzHere) < mzTolHere);
                                // filter by rt
                                binPeaksHere = binPeaksHere.Where(p => Math.Abs(p.retentionTime - identification.ms2RetentionTime) < initialMbrRtWindow);
                                // remove duplicates
                                binPeaksHere = binPeaksHere.Distinct();
                                // filter by isotopic distribution
                                var validIsotopeClusters = FilterPeaksByIsotopicDistribution(binPeaksHere, identification, chargeState, true);

                                if (validIsotopeClusters.Any())
                                {
                                    mbrFeature.isotopeClusters.AddRange(validIsotopeClusters);
                                }
                            }

                            if (mbrFeature.isotopeClusters.Any())
                            {
                                concurrentBagOfMatchedFeatures.Add(mbrFeature);
                            }
                        }
                    }
                );
            }

            features.AddRange(concurrentBagOfMatchedFeatures);
        }

        private void RunErrorChecking(List<ChromatographicPeak> features)
        {
            features.RemoveAll(p => p.isMbrFeature && p.intensity == 0);

            if (!silent)
                Console.WriteLine("Checking errors");
            var featuresWithSamePeak = features.Where(v => v.intensity != 0).GroupBy(p => p.apexPeak.peakWithScan);
            featuresWithSamePeak = featuresWithSamePeak.Where(p => p.Count() > 1);

            // condense duplicate features
            foreach (var duplicateFeature in featuresWithSamePeak)
                duplicateFeature.First().MergeFeatureWith(duplicateFeature, integrate);
            features.RemoveAll(p => p.intensity == -1);

            // check for multiple features per peptide within a time window
            var featuresToMaybeMerge = features.Where(p => p.numIdentificationsByFullSeq == 1 && p.apexPeak != null).GroupBy(p => p.identifyingScans.First().FullSequence).Where(p => p.Count() > 1);
            if (featuresToMaybeMerge.Any())
            {
                foreach (var group in featuresToMaybeMerge)
                {
                    if (idSpecificChargeState)
                    {
                        var group2 = group.ToList().GroupBy(p => p.apexPeak.chargeState).Where(v => v.Count() > 1);

                        foreach (var group3 in group2)
                        {
                            foreach (var feature in group3)
                            {
                                if (feature.intensity != -1)
                                {
                                    var featuresToMerge = group3.Where(p => Math.Abs(p.apexPeak.peakWithScan.retentionTime - feature.apexPeak.peakWithScan.retentionTime) < rtTol && p.intensity != -1);
                                    if (featuresToMerge.Any())
                                        feature.MergeFeatureWith(featuresToMerge, integrate);
                                }
                            }
                        }
                    }
                    else
                    {
                        foreach (var feature in group)
                        {
                            if (feature.intensity != -1)
                            {
                                var featuresToMerge = group.Where(p => Math.Abs(p.apexPeak.peakWithScan.retentionTime - feature.apexPeak.peakWithScan.retentionTime) < rtTol && p.intensity != -1);
                                if (featuresToMerge.Any())
                                    feature.MergeFeatureWith(featuresToMerge, integrate);
                            }
                        }
                    }
                }

                features.RemoveAll(p => p.intensity == -1);
            }

            if (errorCheckAmbiguousMatches)
            {
                // check for multiple peptides per feature
                var scansWithMultipleDifferentIds = features.Where(p => p.numIdentificationsByFullSeq > 1);
                var ambiguousFeatures = scansWithMultipleDifferentIds.Where(p => p.numIdentificationsByBaseSeq > 1);

                // handle ambiguous features
                foreach (var ambiguousFeature in ambiguousFeatures)
                {
                    var msmsIdentsForThisFile = ambiguousFeature.identifyingScans.Where(p => p.fileName.Equals(ambiguousFeature.fileName));

                    if (!msmsIdentsForThisFile.Any())
                    {
                        // mbr matched more than one identification to this peak - cannot resolve
                        ambiguousFeature.intensity = -1;
                    }
                    else
                    {
                        // msms identifications take precident over mbr features
                        ambiguousFeature.identifyingScans = msmsIdentsForThisFile.ToList();
                        ambiguousFeature.ResolveIdentifications();
                    }
                }

                features.RemoveAll(p => p.intensity == -1);
            }
        }

        public List<Peptide> SumFeatures(IEnumerable<ChromatographicPeak> features, bool sumPeptideIntensitiesRegardlessOfModifications)
        {
            List<Peptide> returnList = new List<Peptide>();

            string[] fileNames = new string[filePaths.Length];
            for (int i = 0; i < fileNames.Length; i++)
                fileNames[i] = Path.GetFileNameWithoutExtension(filePaths[i]);
            Peptide.files = fileNames;

            var sequenceToPeaksMatch = new Dictionary<string, List<ChromatographicPeak>>();
            foreach (var feature in features)
            {
                IEnumerable<IGrouping<string, Identification>> seqs;
                if (sumPeptideIntensitiesRegardlessOfModifications)
                    seqs = feature.identifyingScans.GroupBy(p => p.BaseSequence);
                else
                    seqs = feature.identifyingScans.GroupBy(p => p.FullSequence);

                foreach (var seq in seqs)
                {
                    List<ChromatographicPeak> featuresForThisBaseSeq;
                    if (sequenceToPeaksMatch.TryGetValue(seq.Key, out featuresForThisBaseSeq))
                        featuresForThisBaseSeq.Add(feature);
                    else
                        sequenceToPeaksMatch.Add(seq.Key, new List<ChromatographicPeak>() { feature });
                }
            }

            foreach (var sequence in sequenceToPeaksMatch)
            {
                double[] intensitiesByFile = new double[filePaths.Length];
                string[] identificationType = new string[filePaths.Length];
                var thisSeqPerFile = sequence.Value.GroupBy(p => p.fileName);

                for (int i = 0; i < intensitiesByFile.Length; i++)
                {
                    string file = Path.GetFileNameWithoutExtension(filePaths[i]);
                    var featuresForThisSeqAndFile = thisSeqPerFile.Where(p => p.Key.Equals(file)).FirstOrDefault();

                    if (featuresForThisSeqAndFile != null)
                    {
                        if (featuresForThisSeqAndFile.First().isMbrFeature)
                        {
                            identificationType[i] = "MBR";
                            intensitiesByFile[i] = featuresForThisSeqAndFile.Select(p => p.intensity).Max();
                        }
                        else
                        {
                            identificationType[i] = "MSMS";
                            double summedPeakIntensity = featuresForThisSeqAndFile.Sum(p => p.intensity);

                            if (featuresForThisSeqAndFile.Max(p => p.numIdentificationsByBaseSeq) == 1)
                                intensitiesByFile[i] = summedPeakIntensity;
                            else
                            {
                                double ambigPeakIntensity = featuresForThisSeqAndFile.Where(p => p.numIdentificationsByBaseSeq > 1).Sum(v => v.intensity);

                                if ((ambigPeakIntensity / summedPeakIntensity) < 0.3)
                                    intensitiesByFile[i] = featuresForThisSeqAndFile.Select(p => (p.intensity / p.numIdentificationsByBaseSeq)).Sum();
                                else
                                    intensitiesByFile[i] = -1;
                            }
                        }
                    }
                    else
                        identificationType[i] = "";
                }

                if(sumPeptideIntensitiesRegardlessOfModifications)
                    returnList.Add(new Peptide(sequence.Key, string.Join(";", sequence.Value.SelectMany(p => p.identifyingScans).Where(p => p.BaseSequence.Equals(sequence.Key)).First().proteinGroups.Select(p => p.proteinGroupName)), intensitiesByFile, identificationType));
                else
                    returnList.Add(new Peptide(sequence.Key, string.Join(";", sequence.Value.SelectMany(p => p.identifyingScans).Where(p => p.FullSequence.Equals(sequence.Key)).First().proteinGroups.Select(p => p.proteinGroupName)), intensitiesByFile, identificationType));
            }

            return returnList.OrderBy(p => p.Sequence).ToList();
        }

        private IEnumerable<IsotopeCluster> FilterPeaksByIsotopicDistribution(IEnumerable<IndexedMassSpectralPeak> peaks, Identification identification, int chargeState, bool lookForBadIsotope)
        {
            var isotopeClusters = new List<IsotopeCluster>();
            var isotopeMassShifts = baseSequenceToIsotopicDistribution[identification.BaseSequence];

            if (isotopeMassShifts.Count < numIsotopesRequired)
                return isotopeClusters;

            double isotopeMzTol = ((isotopePpmTolerance / 1e6) * identification.monoisotopicMass) / chargeState;

            foreach (var thisPeakWithScan in peaks)
            {
                // calculate theoretical isotopes relative to observed peak
                var theorIsotopeMzs = new double[isotopeMassShifts.Count];
                int isotopicPeakUnitOfPeakZeroIsMono = Convert.ToInt32(thisPeakWithScan.mainPeak.Mz.ToMass(chargeState) - identification.monoisotopicMass);
                var mainpeakMz = thisPeakWithScan.mainPeak.Mz;

                // left of main peak 
                for (int i = 0; i < isotopicPeakUnitOfPeakZeroIsMono; i++)
                    theorIsotopeMzs[i] = mainpeakMz + (isotopeMassShifts[i].Key / chargeState);

                // main peak and right of main peak
                for (int i = isotopicPeakUnitOfPeakZeroIsMono; i < isotopeMassShifts.Count; i++)
                    theorIsotopeMzs[i] = mainpeakMz + (isotopeMassShifts[i].Key / chargeState);
                theorIsotopeMzs = theorIsotopeMzs.OrderBy(p => p).ToArray();

                var lowestMzIsotopePossible = theorIsotopeMzs.First();
                lowestMzIsotopePossible -= isotopeMzTol;
                var highestMzIsotopePossible = theorIsotopeMzs.Last();
                highestMzIsotopePossible += isotopeMzTol;

                // get possible isotope peaks from the peak's scan
                List<MassSpectralPeak> possibleIsotopePeaks = new List<MassSpectralPeak>();

                // go backwards from the peak to find the lowest-mass isotope possible
                int earliestIsotopicPeakIndexPossible = thisPeakWithScan.zeroBasedIndexOfPeakInScan;
                for (int i = earliestIsotopicPeakIndexPossible; i >= 0; i--)
                {
                    if (thisPeakWithScan.scan.MassSpectrum.XArray[i] < lowestMzIsotopePossible)
                        break;
                    earliestIsotopicPeakIndexPossible = i;
                }

                // find the highest-mass isotope possible
                for (int i = earliestIsotopicPeakIndexPossible; i < thisPeakWithScan.scan.MassSpectrum.Size; i++)
                {
                    if (thisPeakWithScan.scan.MassSpectrum.XArray[i] > highestMzIsotopePossible)
                        break;
                    possibleIsotopePeaks.Add(new MassSpectralPeak(thisPeakWithScan.scan.MassSpectrum.XArray[i], thisPeakWithScan.scan.MassSpectrum.YArray[i]));
                }
                
                if (lookForBadIsotope)
                {
                    bool badPeak = false;
                    double prevIsotopePeakMz = (theorIsotopeMzs[0] - (1.003322 / chargeState));

                    for (int i = earliestIsotopicPeakIndexPossible; i > 0; i--)
                    {
                        if (Math.Abs(thisPeakWithScan.scan.MassSpectrum.XArray[i] - prevIsotopePeakMz) < isotopeMzTol)
                            if (thisPeakWithScan.scan.MassSpectrum.YArray[i] / thisPeakWithScan.mainPeak.Intensity > 1.0)
                                badPeak = true;
                        if (thisPeakWithScan.scan.MassSpectrum.XArray[i] < (prevIsotopePeakMz - isotopeMzTol))
                            break;
                    }

                    if (badPeak)
                        continue;
                }

                // isotopic distribution check
                bool isotopeDistributionCheck = false;
                MassSpectralPeak[] isotopePeaks = new MassSpectralPeak[isotopeMassShifts.Count];
                int isotopeIndex = 0;
                double theorIsotopeMz = theorIsotopeMzs[isotopeIndex];
                foreach (var possibleIsotopePeak in possibleIsotopePeaks)
                {
                    if (Math.Abs(possibleIsotopePeak.Mz - theorIsotopeMz) < isotopeMzTol)
                    {
                        // store the good isotope peak
                        isotopePeaks[isotopeIndex] = possibleIsotopePeak;

                        if (isotopeIndex < isotopePeaks.Length - 1)
                        {
                            // look for the next isotope
                            isotopeIndex++;
                            theorIsotopeMz = theorIsotopeMzs[isotopeIndex];
                        }
                    }
                }

                // all isotopes have been looked for - check to see if they've been observed
                if (isotopePeaks.Where(p => p != null).Count() < numIsotopesRequired)
                    continue;

                if (isotopePeaks[0] == null && requireMonoisotopicMass)
                    continue;
                
                if (requireMonoisotopicMass)
                {
                    bool[] requiredIsotopesSeen = new bool[numIsotopesRequired];

                    for (int i = 0; i < numIsotopesRequired; i++)
                    {
                        if (isotopePeaks[i] != null)
                            requiredIsotopesSeen[i] = true;
                        else
                            requiredIsotopesSeen[i] = false;
                    }

                    if (requiredIsotopesSeen.Where(p => p == false).Any())
                        isotopeDistributionCheck = false;
                    else
                        isotopeDistributionCheck = true;
                }
                else
                {
                    // need to look for continuous isotope distribution when monoisotopic mass observation is not required
                    isotopeDistributionCheck = true;
                }

                if (isotopeDistributionCheck)
                {
                    double isotopeClusterIntensity = 0;

                    if (requireMonoisotopicMass)
                    {
                        for (int i = 0; i < isotopePeaks.Length; i++)
                        {
                            if (isotopePeaks[i] != null)
                            {
                                double relIsotopeAbundance = isotopePeaks[i].Intensity / isotopePeaks[0].Intensity;
                                double theorIsotopeAbundance = isotopeMassShifts[i].Value / isotopeMassShifts[0].Value;

                                // impute isotope intensity if it is very different from expected
                                if ((relIsotopeAbundance / theorIsotopeAbundance) < 2.0)
                                    isotopeClusterIntensity += isotopePeaks[i].Intensity;
                                else
                                    isotopeClusterIntensity += theorIsotopeAbundance * isotopePeaks[0].Intensity * 2.0;
                            }
                            else
                                isotopeClusterIntensity += (isotopeMassShifts[i].Value / isotopeMassShifts[0].Value) * isotopePeaks[0].Intensity;
                        }
                    }
                    else
                    {
                        isotopeClusterIntensity = isotopePeaks.Where(p => p != null).Sum(p => p.Intensity);
                    }

                    isotopeClusters.Add(new IsotopeCluster(thisPeakWithScan, chargeState, isotopeClusterIntensity));
                }
            }

            return isotopeClusters;
        }

        private IEnumerable<IndexedMassSpectralPeak> ScanCrawl(IOrderedEnumerable<IndexedMassSpectralPeak> peaksWithScans, int missedScansAllowed, int startingMS1ScanNumber, List<int> ms1ScanNumbers)
        {
            var validPeaksWithScans = new List<IndexedMassSpectralPeak>();

            int lastGoodIndex = ms1ScanNumbers.IndexOf(startingMS1ScanNumber);
            int ms1IndexHere = lastGoodIndex - 1;
            int missedScans = 0;

            foreach (var thisPeakWithScan in peaksWithScans)
            {
                ms1IndexHere = ms1ScanNumbers.IndexOf(thisPeakWithScan.oneBasedScanNumber);
                missedScans += Math.Abs(ms1IndexHere - lastGoodIndex) - 1;

                if (missedScans > missedScansAllowed)
                    break;

                // found a good peak; reset missed scans to 0
                missedScans = 0;
                lastGoodIndex = ms1IndexHere;

                validPeaksWithScans.Add(thisPeakWithScan);
            }

            return validPeaksWithScans;
        }

        private void CutPeak(ChromatographicPeak peak, bool integrate, List<int> ms1ScanNumbers)
        {
            bool cutThisPeak = false;
            IsotopeCluster valleyTimePoint = null;

            if (peak.isotopeClusters.Count() < 5)
                return;

            // find out if we need to split this peak by using the discrimination factor
            var timePointsForApexZ = peak.isotopeClusters.Where(p => p.chargeState == peak.apexPeak.chargeState);
            var leftTimePoints = timePointsForApexZ.Where(p => p.peakWithScan.retentionTime <= peak.apexPeak.peakWithScan.retentionTime).OrderByDescending(v => v.peakWithScan.retentionTime);
            var rightTimePoints = timePointsForApexZ.Where(p => p.peakWithScan.retentionTime >= peak.apexPeak.peakWithScan.retentionTime).OrderBy(v => v.peakWithScan.retentionTime);

            double mind0 = 0.6;

            foreach (var timePoint in rightTimePoints)
            {
                if (valleyTimePoint == null || timePoint.isotopeClusterIntensity < valleyTimePoint.isotopeClusterIntensity)
                    valleyTimePoint = timePoint;

                var timePointsBetweenApexAndThisTimePoint = rightTimePoints.Where(p => p.peakWithScan.retentionTime <= timePoint.peakWithScan.retentionTime).ToList();

                var d0 = (timePoint.isotopeClusterIntensity - valleyTimePoint.isotopeClusterIntensity) / timePoint.isotopeClusterIntensity;
                if (d0 > mind0)
                {
                    var secondValleyTimePoint = timePointsBetweenApexAndThisTimePoint[timePointsBetweenApexAndThisTimePoint.IndexOf(valleyTimePoint) + 1];

                    d0 = (timePoint.isotopeClusterIntensity - secondValleyTimePoint.isotopeClusterIntensity) / timePoint.isotopeClusterIntensity;

                    if (d0 > mind0)
                    {
                        cutThisPeak = true;
                        break;
                    }
                    else
                    {
                        // check for missed scan around valley time point
                        var tpBeforeValleyTimePoint = timePointsBetweenApexAndThisTimePoint[timePointsBetweenApexAndThisTimePoint.IndexOf(valleyTimePoint) - 1];

                        int indexOfTimepointBeforeValleyScan = ms1ScanNumbers.IndexOf(tpBeforeValleyTimePoint.peakWithScan.oneBasedScanNumber);
                        int indexOfValleyScan = ms1ScanNumbers.IndexOf(valleyTimePoint.peakWithScan.oneBasedScanNumber);
                        int indexOfSecondValleyScan = ms1ScanNumbers.IndexOf(secondValleyTimePoint.peakWithScan.oneBasedScanNumber);

                        if (Math.Abs(indexOfValleyScan - indexOfTimepointBeforeValleyScan) > 1)
                        {
                            cutThisPeak = true;
                            break;
                        }
                        else if (Math.Abs(indexOfValleyScan - indexOfSecondValleyScan) > 1)
                        {
                            cutThisPeak = true;
                            break;
                        }
                    }
                }
            }

            if (cutThisPeak == false)
            {
                valleyTimePoint = null;

                foreach (var timePoint in leftTimePoints)
                {
                    if (valleyTimePoint == null || timePoint.isotopeClusterIntensity < valleyTimePoint.isotopeClusterIntensity)
                        valleyTimePoint = timePoint;

                    var timePointsBetweenApexAndThisTimePoint = leftTimePoints.Where(p => p.peakWithScan.retentionTime >= timePoint.peakWithScan.retentionTime).ToList();

                    var d0 = (timePoint.isotopeClusterIntensity - valleyTimePoint.isotopeClusterIntensity) / timePoint.isotopeClusterIntensity;
                    if (d0 > mind0)
                    {
                        var secondValleyTimePoint = timePointsBetweenApexAndThisTimePoint[timePointsBetweenApexAndThisTimePoint.IndexOf(valleyTimePoint) + 1];

                        d0 = (timePoint.isotopeClusterIntensity - secondValleyTimePoint.isotopeClusterIntensity) / timePoint.isotopeClusterIntensity;

                        if (d0 > mind0)
                        {
                            cutThisPeak = true;
                            break;
                        }
                        else
                        {
                            // check for missed scan around valley time point
                            var tpBeforeValleyTimePoint = timePointsBetweenApexAndThisTimePoint[timePointsBetweenApexAndThisTimePoint.IndexOf(valleyTimePoint) - 1];

                            int indexOfTimepointBeforeValleyScan = ms1ScanNumbers.IndexOf(tpBeforeValleyTimePoint.peakWithScan.oneBasedScanNumber);
                            int indexOfValleyScan = ms1ScanNumbers.IndexOf(valleyTimePoint.peakWithScan.oneBasedScanNumber);
                            int indexOfSecondValleyScan = ms1ScanNumbers.IndexOf(secondValleyTimePoint.peakWithScan.oneBasedScanNumber);

                            if (Math.Abs(indexOfValleyScan - indexOfTimepointBeforeValleyScan) > 1)
                            {
                                cutThisPeak = true;
                                break;
                            }
                            else if (Math.Abs(indexOfValleyScan - indexOfSecondValleyScan) > 1)
                            {
                                cutThisPeak = true;
                                break;
                            }
                        }
                    }
                }
            }

            // cut
            if (cutThisPeak)
            {
                var splitLeft = peak.isotopeClusters.Where(p => p.peakWithScan.retentionTime <= valleyTimePoint.peakWithScan.retentionTime).ToList();
                var splitRight = peak.isotopeClusters.Where(p => p.peakWithScan.retentionTime >= valleyTimePoint.peakWithScan.retentionTime).ToList();

                if (peak.identifyingScans.First().ms2RetentionTime > splitLeft.Max(p => p.peakWithScan.retentionTime))
                    foreach (var timePoint in splitLeft)
                        peak.isotopeClusters.Remove(timePoint);
                else
                    foreach (var timePoint in splitRight)
                        peak.isotopeClusters.Remove(timePoint);

                // recalculate intensity for the peak
                peak.CalculateIntensityForThisFeature(integrate);
                peak.splitRT = valleyTimePoint.peakWithScan.retentionTime;

                // recursively cut
                CutPeak(peak, integrate, ms1ScanNumbers);
            }
        }
    }
}