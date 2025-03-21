using System.Text;
using Omics.Modifications;
using Transcriptomics.Digestion;

namespace Transcriptomics
{
    public static class ClassExtensions
    {
        public static T CreateNew<T>(this T target, string? sequence = null, IDictionary<int, List<Modification>>? modifications = null,
        bool? isDecoy = null)
            where T : INucleicAcid
        {
            // set new object parameters where not null
            object? returnObj = null;
            string newSequence = sequence ?? target.BaseSequence;
            IDictionary<int, List<Modification>> newModifications = modifications ?? target.OneBasedPossibleLocalizedModifications;


            switch (target)
            {
                case RNA rna:
                {
                    bool newIsDecoy = isDecoy ?? rna.IsDecoy;
                    returnObj = new RNA(newSequence, rna.Name, rna.Accession, rna.Organism, rna.DatabaseFilePath,
                        rna.FivePrimeTerminus, rna.ThreePrimeTerminus, newModifications, rna.IsContaminant, newIsDecoy, rna.GeneNames.ToList(), rna.AdditionalDatabaseFields);
                    break;
                }
                case OligoWithSetMods oligo:
                {
                    var oldParent = oligo.Parent as RNA ?? throw new NullReferenceException();
                    bool newIsDecoy = isDecoy ?? oldParent.IsDecoy;
                    var newParent = new RNA(
                        newSequence,
                        oldParent.Name,
                        oldParent.Accession,
                        oldParent.Organism,
                        oldParent.DatabaseFilePath,
                        oldParent.FivePrimeTerminus,
                        oldParent.ThreePrimeTerminus,
                        newModifications,
                        oldParent.IsContaminant,
                        newIsDecoy,
                        oldParent.GeneNames.ToList(),
                        oldParent.AdditionalDatabaseFields);

                    returnObj = new OligoWithSetMods(
                        newParent,
                        (oligo.DigestionParams as RnaDigestionParams)!,
                        oligo.OneBasedStartResidue,
                        oligo.OneBasedEndResidue,
                        oligo.MissedCleavages,
                        oligo.CleavageSpecificityForFdrCategory,
                        newModifications.ToDictionary(p => p.Key, p => p.Value.First()),
                        oligo.NumFixedMods,
                        oligo.FivePrimeTerminus,
                        oligo.ThreePrimeTerminus);
                    break;
                }
                default:
                    throw new ArgumentException("INucleicAcid type not yet implemented");
            }

            return (T)returnObj ?? throw new NullReferenceException("Error creating new INucleicAcid");
        }

        /// <summary>
        /// Transcribes a DNA sequence into an RNA sequence
        /// </summary>
        /// <param name="dna">The input dna sequence</param>
        /// <param name="isCodingStrand">True if the input sequence is the coding strand, False if the input sequence is the template strand</param>
        /// <returns></returns>
        public static string Transcribe(this string dna, bool isCodingStrand = true)
        {
            var sb = new StringBuilder();
            foreach (var residue in dna)
            {
                if (isCodingStrand)
                {
                    sb.Append(residue == 'T' ? 'U' : residue);
                }
                else
                {
                    switch (residue)
                    {
                        case 'A':
                            sb.Append('U');
                            break;
                        case 'T':
                            sb.Append('A');
                            break;
                        case 'C':
                            sb.Append('G');
                            break;
                        case 'G':
                            sb.Append('C');
                            break;
                        default:
                            sb.Append(residue);
                            break;
                    }
                }
            }
            return sb.ToString();
        }
    }
}
