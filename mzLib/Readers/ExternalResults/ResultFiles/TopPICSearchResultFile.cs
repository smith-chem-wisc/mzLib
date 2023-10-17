using CsvHelper.Configuration;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Reflection.Metadata.Ecma335;
using System.Text;
using System.Threading.Tasks;
using Easy.Common.Extensions;
using MassSpectrometry;
using ThermoFisher.CommonCore.Data;
using CsvHelper;
using CsvHelper.TypeConversion;
using Microsoft.SqlServer.Server;
using Readers.ExternalResults.IndividualResultRecords;

namespace Readers
{
    public class ToppicSearchResultFile : ResultFile<ToppicPrsm>
    {
        public override SupportedFileType FileType
        {
            get
            {
                if (FilePath.EndsWith(SupportedFileType.ToppicProteoform.GetFileExtension()))
                    return SupportedFileType.ToppicProteoform;
                if (FilePath.EndsWith(SupportedFileType.ToppicPrsm.GetFileExtension()))
                    return SupportedFileType.ToppicPrsm;
                if (FilePath.EndsWith(SupportedFileType.ToppicProteoformSingle.GetFileExtension()))
                    return SupportedFileType.ToppicProteoformSingle;
                if (FilePath.EndsWith(SupportedFileType.ToppicPrsmSingle.GetFileExtension()))
                    return SupportedFileType.ToppicPrsmSingle;
                throw new FileLoadException("File type not supported");
            }
        }
        public override Software Software { get; set; }

        #region Search Summary Properties

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
                    bool isAlternativeID = false;
                    var toppicDefaultConfig = ToppicPrsm.CsvConfiguration;
                    var csvConfig = new CsvConfiguration(CultureInfo.InvariantCulture)
                    {
                        Delimiter = toppicDefaultConfig.Delimiter,
                        Encoding = toppicDefaultConfig.Encoding,
                        HasHeaderRecord = toppicDefaultConfig.HasHeaderRecord,
                        ReadingExceptionOccurred = context =>
                        {
                            if (context.Exception is not TypeConverterException) throw new IOException("Error reading Toppic results file", context.Exception);

                            isAlternativeID = true;
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
                        if (isAlternativeID)
                        {
                            results.Last().AlternativeIdentifications.AddRange(alternativeIDs
                                .Select(p => p.Split('\t')
                                    .Where(str => !string.IsNullOrEmpty(str))
                                    .ToArray())
                                .Select(p => (int.Parse(p[1]), p[2], p[3], int.Parse(p[4]), int.Parse(p[5])))
                                .ToList());
                            isAlternativeID = false;
                            alternativeIDs.Clear();
                        }
                        else
                        {
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

            using var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), ToppicPrsm.CsvConfiguration);

            csv.WriteHeader<ToppicPrsm>();
            foreach (var result in Results)
            {
                csv.NextRecord();
                csv.WriteRecord(result);
            }
        }

    
    }
}
