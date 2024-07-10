using CsvHelper.Configuration.Attributes;
using CsvHelper.Configuration;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers
{
    public class MsFraggerExperiment
    {
        public static CsvConfiguration CsvConfiguration = new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Delimiter = "\t",
            HasHeaderRecord = true,
            IgnoreBlankLines = true,
            TrimOptions = TrimOptions.Trim,
            BadDataFound = null,
            MissingFieldFound = null,
            HeaderValidated = null,
        };

        [Name("file")] public string FullFilePathWithExtension { get; set; }
        [Name("sample")][Optional] public string Sample { get; set; }

        [Name("sample_name")][Optional] public string SampleName { get; set; }

        [Name("condition")][Optional] public string Condition { get; set; }
        [Name("replicate")][Optional] public int Replicate { get; set; }   

        public MsFraggerExperiment()
        {
            FullFilePathWithExtension = string.Empty;
            Sample = string.Empty;
            SampleName = string.Empty;
            Condition = string.Empty;
            Replicate = 0;
        }
    }
}
