using CsvHelper.Configuration.Attributes;
using CsvHelper.Configuration;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers
{
    public class MsFraggerProtein
    {
        public static CsvConfiguration CsvConfiguration => new CsvConfiguration(System.Globalization.CultureInfo.InvariantCulture)
        {
            Delimiter = "\t",
            HasHeaderRecord = true,
            IgnoreBlankLines = true,
            TrimOptions = TrimOptions.Trim,
            BadDataFound = null,
            MissingFieldFound = null,
        };

        [Name("Protein")]
        public string Protein { get; set; }

        [Name("Protein ID")]
        public string Accesion { get; set; }

        [Name("Entry Name")]
        public string AccessionOrganism { get; set; }

        [Name("Gene")]
        public string Gene { get; set; }

        [Name("Length", "Protein Length")]
        public int Length { get; set; }

        [Name("Organism")]
        public string Organism { get; set; }

        [Name("Protein Description", "Description")]
        public string Description { get; set; }

        [Name("Protein Existence")]
        public string ProteinExistence { get; set; }

        [Name("Coverage")]
        [Optional]
        public double Coverage { get; set; }

        [Name("Protein Probability")]
        public double ProteinProbability { get; set; }

        [Name("Top Peptide Probability")]
        public double TopPeptideProbability { get; set; }

        [Name("Total Peptides", "Combined Total Peptides")]
        public int TotalPeptides { get; set; }

        [Name("Unique Peptides")]
        [Optional]
        public int UniquePeptides { get; set; }

        [Name("Razor Peptides")]
        [Optional]
        public int RazorPeptides { get; set; }

        [Name("Total Spectral Count", "Combined Total Spectral Count")]
        public int TotalSpectralCount { get; set; }

        [Name("Unique Spectral Count", "Combined Unique Spectral Count")]
        public int UniqueSpectralCount { get; set; }

        [Name("Razor Spectral Count")]
        [Optional]
        public int RazorSpectralCount { get; set; }

        [Name("Total Intensity")]
        [Optional]
        public double TotalIntensity { get; set; }

        [Name("Unique Intensity")]
        [Optional]
        public double UniqueIntensity { get; set; }

        [Name("Razor Intensity")]
        [Optional]
        public double RazorIntensity { get; set; }

        [Name("Razor Assigned Modifications")]
        [TypeConverter(typeof(CommaDelimitedToStringArrayTypeConverter))]
        [Optional]
        public string[] RazorAssignedModifications { get; set; }

        [Name("Razor Observed Modifications")]
        [Optional]
        [TypeConverter(typeof(CommaDelimitedToStringArrayTypeConverter))]
        public string[] RazorObservedModifications { get; set; }

        [Name("Indistinguishable Proteins")]
        [TypeConverter(typeof(CommaDelimitedToStringArrayTypeConverter))]
        public string[] IndistinguishableProteins { get; set; }

        public MsFraggerProtein()
        {
            RazorAssignedModifications = new string[0];
            RazorObservedModifications = new string[0];
            IndistinguishableProteins = new string[0];
        }
    }
}
