using CsvHelper.Configuration;
using System.Globalization;
using System.Text;
using Easy.Common.Extensions;
using MassSpectrometry;
using CsvHelper;
using CsvHelper.TypeConversion;
using MzLibUtil;

namespace Readers
{
    /// <summary>
    /// Concrete Product for reading and representing a proteoform or psm search results file from TopPIC
    /// For supported versions and software this file type can come from see
    ///     Readers.ExternalResources.SupportedVersions.txt
    /// </summary>
    public class ToppicSearchResultFile : ResultFile<ToppicPrsm>
    {
        private SupportedFileType _fileType;
        public override SupportedFileType FileType
        {
            get
            {
                if (!_fileType.IsDefault()) return _fileType;

                if (FilePath.EndsWith(SupportedFileType.ToppicProteoform.GetFileExtension()))
                    _fileType = SupportedFileType.ToppicProteoform;
                else if (FilePath.EndsWith(SupportedFileType.ToppicPrsm.GetFileExtension()))
                    _fileType = SupportedFileType.ToppicPrsm;
                else if (FilePath.EndsWith(SupportedFileType.ToppicProteoformSingle.GetFileExtension()))
                    _fileType = SupportedFileType.ToppicProteoformSingle;
                else if (FilePath.EndsWith(SupportedFileType.ToppicPrsmSingle.GetFileExtension()))
                    _fileType = SupportedFileType.ToppicPrsmSingle;
                else throw new MzLibException("Cannot parse result file type from file path");

                return _fileType;
            }
        }
        public override Software Software { get; set; }

        #region Search Summary Parameters

        public string ProteinDatabasePath { get; private set; }
        public string SpectrumFilePath { get; private set; }
        public int NumberOfCombinedSpectra { get; private set; }
        public DissociationType FragmentationMethod { get; private set; }
        public string SearchType { get; private set; }
        public List<string> FixedModifications { get; private set; }
        public List<string> AllowedNTerminalForms { get; private set; }
        public int NumberOfMaxUnexpectedModifications { get; private set; }
        public double MaximumMassShift { get; private set; }
        public double MinimumMassShift { get; private set; }
        public string SpectrumLevelCutOffType { get; private set; }
        public double SpectrumLevelCutOffValue { get; private set; }
        public string ProteoformLevelCutOffType { get; private set; }
        public double ProteoformLevelCutOffValue { get; private set; }
        public double PrecursorErrorTolerance { get; private set; }
        public double PrsmClusterErrorTolerance { get; private set; }
        public bool UseToppicFeatureFile { get; private set; }
        public string EValueComputation { get; private set; }
        public bool LocalizationWithMIScore { get; private set; }
        public int ThreadNumber { get; private set; }
        public string ExecutableFileDirectory { get; private set; }
        public DateTime StartTime { get; private set; }
        public DateTime EndTime { get; private set; }
        public string Version { get; private set; }

        #endregion

        public ToppicSearchResultFile(string filePath) : base(filePath, Software.Toppic)
        {
            FixedModifications = new List<string>();
            AllowedNTerminalForms = new List<string>();
        }

        /// <summary>
        /// Constructor used to initialize from the factory method
        /// </summary>
        public ToppicSearchResultFile() : base()
        {
            FixedModifications = new List<string>();
            AllowedNTerminalForms = new List<string>();
        }


        public override void LoadResults()
        {
            using var reader = new StreamReader(FilePath);

            // Pull out parameter values
            bool isInsideParametersSection = false;
            bool isInsideFixedModificationsSection = false;
            StringBuilder dataBetweenParameters = new StringBuilder();

            // Read the file line by line
            while (reader.ReadLine() is { } line)
            {
                if (line.Contains("******* Parameters *******"))
                {
                    // Found the start or end of the Parameters section
                    if (isInsideParametersSection)
                    {
                        // We are done, with parameters
                        isInsideParametersSection = false;
                        continue;
                    }
                    else
                    {
                        // We are entering the Parameters section
                        isInsideParametersSection = true;
                        continue; // Skip this line
                    }
                }
                if (line.Equals(""))
                {
                    // everything below will be the actual Toppic results
                    // read them in with CsvHelper
                    var alternativeIDs = new List<string>();
                    bool isAlternativeId = false;
                    var toppicDefaultConfig = ToppicPrsm.CsvConfiguration;
                    var csvConfig = new CsvConfiguration(CultureInfo.InvariantCulture)
                    {
                        Delimiter = toppicDefaultConfig.Delimiter,
                        Encoding = toppicDefaultConfig.Encoding,
                        HasHeaderRecord = toppicDefaultConfig.HasHeaderRecord,
                        ReadingExceptionOccurred = context =>
                        {
                            if (context.Exception is not TypeConverterException) throw new IOException("Error reading Toppic results file", context.Exception);

                            isAlternativeId = true;
                            alternativeIDs.Add(context.Exception.Context.Parser.RawRecord);
                            return false;
                        },
                    };

                    var results = new List<ToppicPrsm>();
                    using var csv = new CsvReader(reader, csvConfig);

                    csv.Read();
                    csv.ReadHeader();
                    while (csv.Read())
                    {
                        var record = csv.GetRecord<ToppicPrsm>();
                        if (isAlternativeId)
                        {
                            results.Last().AlternativeIdentifications.AddRange(alternativeIDs
                                .Select(p => p.Split('\t')
                                    .Where(str => !string.IsNullOrEmpty(str))
                                    .ToArray())
                                .Select(p => new AlternativeToppicId(int.Parse(p[1]), p[2], p[3], int.Parse(p[4]), int.Parse(p[5])))
                                .ToList());
                            isAlternativeId = false;
                            alternativeIDs.Clear();
                        }
                        else
                        {
                            if (record != null)
                                results.Add(record);
                        }
                    }

                    Results = results;

                    // dispose of used readers and exit the loop
                    csv.Dispose();
                    reader.Dispose();
                    break;
                }
                if (isInsideParametersSection)
                {
                    if (line.Contains("Fixed modifications BEGIN"))
                    {
                        isInsideFixedModificationsSection = true;
                        continue; // Skip this line
                    }

                    if (isInsideFixedModificationsSection)
                    {
                        if (line.Contains("Fixed modifications END"))
                        {
                            isInsideFixedModificationsSection = false;
                            continue; // Skip this line
                        }
                        else
                        {
                            FixedModifications.Add(string.Join(" ", line.SplitAndTrim('\t', ' ').Where(p => p.IsNotNullOrEmptyOrWhiteSpace())));
                            continue; // Skip this line
                        }
                    }

                    // Append the line to the data between the Parameters sections
                    dataBetweenParameters.AppendLine(line);
                }
            }

            // Display or process the data between the Parameters sections
            var parameterResults = dataBetweenParameters.ToString();

            // Parse the parameter values
            var parameterLines = parameterResults.Split('\n')
                .Where(p => !string.IsNullOrWhiteSpace(p) && !p.Contains("****"))
                .Select(p => p.Split('\t'))
                .ToDictionary(p => p[0].Trim().Replace(":",""), p => p[1].Replace("\r",""));

            string dateTimeFormat = "ddd MMM dd HH:mm:ss yyyy";
            foreach (var parameter in parameterLines)
            {
                switch (parameter.Key)
                {
                    case "Protein database file":
                        ProteinDatabasePath = parameter.Value;
                        break;
                    case "Spectrum file":
                        SpectrumFilePath = parameter.Value;
                        break;
                    case "Number of combined spectra":
                        NumberOfCombinedSpectra = int.Parse(parameter.Value);
                        break;
                    case "Fragmentation method":
                        FragmentationMethod = Enum.Parse<DissociationType>(parameter.Value);
                        break;
                    case "Search type":
                        SearchType = parameter.Value;
                        break;
                    case "Allowed N-terminal forms":
                        AllowedNTerminalForms = parameter.Value.Split(',').ToList();
                        break;
                    case "Maximum number of unexpected modifications":
                        NumberOfMaxUnexpectedModifications = int.Parse(parameter.Value);
                        break;
                    case "Maximum mass shift of modifications":
                        MaximumMassShift = double.Parse(parameter.Value.Split(" ")[0]);
                        break;
                    case "Minimum mass shift of modifications":
                        MinimumMassShift = double.Parse(parameter.Value.Split(" ")[0]);
                        break;
                    case "Spectrum-level cutoff type":
                        SpectrumLevelCutOffType = parameter.Value;
                        break;
                    case "Spectrum-level cutoff value":
                        SpectrumLevelCutOffValue = double.Parse(parameter.Value);
                        break;
                    case "Proteoform-level cutoff type":
                        ProteoformLevelCutOffType = parameter.Value;
                        break;
                    case "Proteoform-level cutoff value":
                        ProteoformLevelCutOffValue = double.Parse(parameter.Value);
                        break;
                    case "Error tolerance for matching masses":
                        PrecursorErrorTolerance = double.Parse(parameter.Value.Split(" ")[0]);
                        break;
                    case "Error tolerance for identifying PrSM clusters":
                        PrsmClusterErrorTolerance = double.Parse(parameter.Value.Split(" ")[0]);
                        break;
                    case "Use TopFD feature file":
                        UseToppicFeatureFile = bool.Parse(parameter.Value);
                        break;
                    case "E-value computation":
                        EValueComputation = parameter.Value;
                        break;
                    case "Localization with MIScore":
                        LocalizationWithMIScore = bool.Parse(parameter.Value);
                        break;
                    case "Thread number":
                        ThreadNumber = int.Parse(parameter.Value);
                        break;
                    case "Executable file directory":
                        ExecutableFileDirectory = parameter.Value;
                        break;
                    case "Start time":
                        if (DateTime.TryParseExact(parameter.Value, dateTimeFormat,
                                System.Globalization.CultureInfo.InvariantCulture,
                                System.Globalization.DateTimeStyles.None, out DateTime result))
                        {
                            StartTime = result;
                        }
                        break;
                    case "End time":
                        if (DateTime.TryParseExact(parameter.Value, dateTimeFormat,
                                System.Globalization.CultureInfo.InvariantCulture,
                                System.Globalization.DateTimeStyles.None, out result))
                        {
                            EndTime = result;
                        }
                        break;
                    case "Version":
                        Version = parameter.Value;
                        break;
                    default:
                        throw new NotImplementedException();
                }
            }
        }

        public override void WriteResults(string outputPath)
        {
            if (!CanRead(outputPath))
                outputPath += FileType.GetFileExtension();

            using var sw = new StreamWriter(File.Create(outputPath));
            using var csv = new CsvWriter(sw, ToppicPrsm.CsvConfiguration);

            // write header
            sw.WriteLine("********************** Parameters **********************");
            sw.WriteLine($"{"Protein database file:",-46}\t{ProteinDatabasePath}");
            sw.WriteLine($"{"Spectrum file:",-46}\t{SpectrumFilePath}");
            sw.WriteLine($"{"Number of combined spectra:",-46}\t{NumberOfCombinedSpectra}");
            sw.WriteLine($"{"Fragmentation method:",-46}\t{FragmentationMethod}");
            sw.WriteLine($"{"Search type:",-46}\t{SearchType}");
            sw.WriteLine($"{"Fixed modifications BEGIN",-46}");
            foreach (var fixedMod in FixedModifications)
            {
                var splits = fixedMod.Split(' ');
                sw.WriteLine($"{splits[0],-46}\t{splits[1]}\t{splits[2]}");
            }
            sw.WriteLine($"{"Fixed modifications END",-46}");
            sw.WriteLine($"{"Allowed N-terminal forms:",-46}\t{string.Join(",", AllowedNTerminalForms)}");
            sw.WriteLine($"{"Maximum number of unexpected modifications:",-46}\t{NumberOfMaxUnexpectedModifications}");
            sw.WriteLine($"{"Maximum mass shift of modifications:",-46}\t{MaximumMassShift} Da");
            sw.WriteLine($"{"Minimum mass shift of modifications:",-46}\t{MinimumMassShift} Da");
            sw.WriteLine($"{"Spectrum-level cutoff type:",-46}\t{SpectrumLevelCutOffType}");
            sw.WriteLine($"{"Spectrum-level cutoff value:",-46}\t{SpectrumLevelCutOffValue}");
            sw.WriteLine($"{"Proteoform-level cutoff type:",-46}\t{ProteoformLevelCutOffType}");
            sw.WriteLine($"{"Proteoform-level cutoff value:",-46}\t{ProteoformLevelCutOffValue}");
            sw.WriteLine($"{"Error tolerance for matching masses:",-46}\t{PrecursorErrorTolerance} ppm");
            sw.WriteLine($"{"Error tolerance for identifying PrSM clusters:",-46}\t{PrsmClusterErrorTolerance} Da");
            sw.WriteLine($"{"Use TopFD feature file:",-46}\t{UseToppicFeatureFile}");
            sw.WriteLine($"{"E-value computation:",-46}\t{EValueComputation}");
            sw.WriteLine($"{"Localization with MIScore:",-46}\t{LocalizationWithMIScore}");
            sw.WriteLine($"{"Thread number:",-46}\t{ThreadNumber}");
            sw.WriteLine($"{"Executable file directory:",-46}\t{ExecutableFileDirectory}");
            sw.WriteLine($"{"Start time:",-46}\t{StartTime:ddd MMM dd HH:mm:ss yyyy}");
            sw.WriteLine($"{"End time:",-46}\t{EndTime:ddd MMM dd HH:mm:ss yyyy}");
            sw.WriteLine($"{"Version:",-46}\t{Version}");
            sw.WriteLine("********************** Parameters **********************");
            sw.WriteLine("");

            csv.WriteHeader<ToppicPrsm>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
                foreach (var alternativeId in result.AlternativeIdentifications)
                {
                    csv.NextRecord();
                    sw.Write($"{result.FilePath}\t{alternativeId}");
                }
            }
        }
    }
}
