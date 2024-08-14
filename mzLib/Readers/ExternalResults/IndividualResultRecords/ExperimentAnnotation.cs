using CsvHelper.Configuration;
using CsvHelper.Configuration.Attributes;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers
{
    public class ExperimentAnnotation
    {
        public static CsvConfiguration CsvConfiguration = new CsvConfiguration(CultureInfo.InvariantCulture)
        {
            Delimiter = "\t",
            HasHeaderRecord = true,
            IgnoreBlankLines = true,
            TrimOptions = TrimOptions.Trim,
            BadDataFound = null,
        };

        #region experiment_annotation Fields

        [Name("file")]
        public string File { get; set; }

        [Name("sample")]
        public string Sample { get; set; }

        [Name("sample_name")]
        public string SampleName { get; set; }

        [Name("condition")]
        public string Condition { get; set; }

        [Name("replicate")]
        public string Replicate { get; set; }

        #endregion
    }
}
