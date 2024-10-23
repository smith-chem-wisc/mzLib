using Chemistry;
using Omics.Modifications;

namespace Transcriptomics
{
    public class RNA : NucleicAcid
    {
        /// <summary>
        /// For constructing RNA from a string
        /// </summary>
        /// <param name="sequence"></param>
        /// <param name="fivePrimeTerm"></param>
        /// <param name="threePrimeTerm"></param>
        /// <param name="oneBasedPossibleLocalizedModifications"></param>
        public RNA(string sequence, IHasChemicalFormula? fivePrimeTerm = null, IHasChemicalFormula? threePrimeTerm = null,
            IDictionary<int, List<Modification>>? oneBasedPossibleLocalizedModifications = null)
            : base(sequence, fivePrimeTerm, threePrimeTerm, oneBasedPossibleLocalizedModifications)
        {
        }

        /// <summary>
        /// For use with RNA loaded from a database
        /// </summary>
        /// <param name="sequence"></param>
        /// <param name="name"></param>
        /// <param name="identifier"></param>
        /// <param name="organism"></param>
        /// <param name="databaseFilePath"></param>
        /// <param name="fivePrimeTerminus"></param>
        /// <param name="threePrimeTerminus"></param>
        /// <param name="oneBasedPossibleModifications"></param>
        /// <param name="isContaminant"></param>
        /// <param name="isDecoy"></param>
        /// <param name="geneNames"></param>
        /// <param name="databaseAdditionalFields"></param>
        public RNA(string sequence, string name, string identifier, string organism, string databaseFilePath,
            IHasChemicalFormula? fivePrimeTerminus = null, IHasChemicalFormula? threePrimeTerminus = null,
            IDictionary<int, List<Modification>>? oneBasedPossibleModifications = null,
            bool isContaminant = false, bool isDecoy = false, List<Tuple<string, string>> geneNames = null,
            Dictionary<string, string>? databaseAdditionalFields = null)
            : base(sequence, name, identifier, organism, databaseFilePath, fivePrimeTerminus, threePrimeTerminus,
                oneBasedPossibleModifications, isContaminant, isDecoy, geneNames, databaseAdditionalFields)
        {

        }
    }
}
