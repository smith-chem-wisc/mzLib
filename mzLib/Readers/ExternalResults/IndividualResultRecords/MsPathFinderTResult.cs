using CsvHelper.Configuration.Attributes;
using CsvHelper.Configuration;
using Chemistry;
using Omics.Modifications;
using Omics;

namespace Readers
{
    public class MsPathFinderTResult : ISpectralMatch
    {
        public static CsvConfiguration CsvConfiguration { get; } = new CsvConfiguration(System.Globalization.CultureInfo.InvariantCulture)
        {
            Delimiter = "\t",
            HasHeaderRecord = true,
            IgnoreBlankLines = true,
            TrimOptions = CsvHelper.Configuration.TrimOptions.InsideQuotes,
            BadDataFound = null,
        };

        [Name("Scan")]
        public int OneBasedScanNumber { get; set; }

        [Name("Pre")]
        public char PreviousResidue { get; set; }

        [Name("Sequence")]
        public string BaseSequence { get; set; }

        [Name("Post")]
        public char NextResidue { get; set; }

        [Name("Modifications")]
        public string Modifications { get; set; }

        [Name("Composition")]
        [TypeConverter(typeof(MsPathFinderTCompositionToChemicalFormulaConverter))]
        public ChemicalFormula ChemicalFormula { get; set; }

        [Name("ProteinName")]
        public string ProteinName { get; set; }

        [Name("ProteinDesc")]
        public string ProteinDescription { get; set; }

        [Name("ProteinLength")]
        public int Length { get; set; }

        [Name("Start")]
        public int OneBasedStartResidue { get; set; }

        [Name("End")]
        public int OneBasedEndResidue { get; set; }

        [Name("Charge")]
        public int Charge { get; set; }

        [Name("MostAbundantIsotopeMz")]
        public double MostAbundantIsotopeMz { get; set; }

        [Name("Mass")]
        public double MonoisotopicMass { get; set; }

        [Name("Ms1Features")]
        public int Ms1Features { get; set; }

        [Name("#MatchedFragments")]
        public int NumberOfMatchedFragments { get; set; }

        [Name("Probability")]
        public double Probability { get; set; }

        [Name("SpecEValue")]
        public double SpecEValue { get; set; }

        [Name("EValue")]
        public double EValue { get; set; }

        [Name("QValue")]
        [Optional]
        public double QValue { get; set; }

        [Name("PepQValue")]
        [Optional]
        public double PepQValue { get; set; }

        #region InterpretedFields

        [Ignore] private string? _accession = null;
        [Ignore] public string Accession => _accession ??= ProteinName.Split('|')[1].Trim();

        [Ignore] private bool? _isDecoy = null;
        [Ignore] public bool IsDecoy => _isDecoy ??= ProteinName.StartsWith("XXX");
        [Optional] public string FileNameWithoutExtension { get; set; }

        [Ignore] private Dictionary<int, Modification>? _allModsOneIsNterminus;
        [Ignore] public Dictionary<int, Modification> AllModsOneIsNterminus => _allModsOneIsNterminus ??= ParseModifications();

        [Ignore] private string? _fullSequence;
        [Ignore] public string FullSequence => _fullSequence ??= IBioPolymerWithSetMods.DetermineFullSequence(BaseSequence, AllModsOneIsNterminus);

        private Dictionary<int, Modification> ParseModifications()
        {
            var mods = new Dictionary<int, Modification>();
            if (string.IsNullOrEmpty(Modifications))
                return mods;

            var modStrings = Modifications.Split(',');
            foreach (var modString in modStrings)
            {
                var modSplits = modString.Split(' ');
                var name = modSplits[0];
                var location = int.Parse(modSplits[1]);

                var modifiedResidue = location == 0 ? 'X' : BaseSequence[location - 1];
                var mmMod = ModificationConverter.GetClosestMod(name, modifiedResidue);

                mods.Add(location+1, mmMod);
            }
            return mods;
        }

        #endregion
    }
}
