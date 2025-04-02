using CsvHelper.Configuration.Attributes;
using CsvHelper.Configuration;
using Chemistry;
using System.Text;

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
        [TypeConverter(typeof(MsPathFinderTPsmStringToModificationsArrayConverter))]
        public MsPathFinderTModification[] Modifications { get; set; }

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

        [Ignore] private string? _fullSequence;
        [Ignore]
        public string FullSequence
        {
            get
            {
                if (_fullSequence != null)
                    return _fullSequence;
                if (!Modifications.Any())
                    return _fullSequence = BaseSequence;

                var sb = new StringBuilder();

                if (Modifications.Any(p => p.OneBasedLocalization == 0))
                {
                    ILocalizedModification? modToAdd = Modifications.FirstOrDefault(p => p.OneBasedLocalization == 0);
                    if (modToAdd is not null)
                        sb.Append(modToAdd.GetMetaMorpheusFullSequenceString());
                }
                for (int i = 0; i < BaseSequence.Length; i++)
                {
                    var residue = BaseSequence[i];
                    sb.Append(residue);

                    ILocalizedModification? potentialMod = Modifications.FirstOrDefault(p => p.OneBasedLocalization == i + 1);
                    if (potentialMod is null) continue;

                    var mmMod = potentialMod.GetMetaMorpheusFullSequenceString();
                    sb.Append(mmMod);
                }

                return _fullSequence = sb.ToString();
            }
        }

        #endregion
    }
}
