using System.Collections.Generic;

namespace Proteomics
{
    public class Protein
    {

        #region Public Constructors

        public Protein(string sequence, string accession, string name, string full_name, bool isDecoy, bool isContaminant)
        {
            BaseSequence = sequence;
            Accession = accession;
            Name = name;
            FullName = full_name;
            IsDecoy = isDecoy;
            IsContaminant = isContaminant;
        }

        public Protein(string sequence, string accession, IDictionary<int, List<Modification>> oneBasedModifications, int?[] oneBasedBeginPositionsForProteolysisProducts, int?[] oneBasedEndPositionsForProteolysisProducts, string[] oneBasedProteolysisProductsTypes, string name, string full_name, bool isDecoy, bool isContaminant, List<GoTerm> goTerms)
        : this(sequence, accession, name, full_name, isDecoy, isContaminant)
        {
            var proteolysisProducts = new List<ProteolysisProduct>();
            if (oneBasedProteolysisProductsTypes != null
                && oneBasedEndPositionsForProteolysisProducts != null
                && oneBasedEndPositionsForProteolysisProducts != null
                && oneBasedProteolysisProductsTypes.Length == oneBasedBeginPositionsForProteolysisProducts.Length
                && oneBasedProteolysisProductsTypes.Length == oneBasedEndPositionsForProteolysisProducts.Length)
                for (int i = 0; i < oneBasedProteolysisProductsTypes.Length; i++)
                    proteolysisProducts.Add(new ProteolysisProduct(oneBasedBeginPositionsForProteolysisProducts[i],
                                                                   oneBasedEndPositionsForProteolysisProducts[i],
                                                                   oneBasedProteolysisProductsTypes[i]));
            ProteolysisProducts = proteolysisProducts;
            GoTerms = goTerms;
            OneBasedPossibleLocalizedModifications = oneBasedModifications;
        }

        #endregion Public Constructors

        #region Public Properties

        public IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; private set; }
        public string Accession { get; private set; }
        public string BaseSequence { get; private set; }
        public bool IsDecoy { get; private set; }
        public IEnumerable<ProteolysisProduct> ProteolysisProducts { get; private set; }
        public IEnumerable<GoTerm> GoTerms { get; private set; }

        public int Length
        {
            get
            {
                return BaseSequence.Length;
            }
        }

        public string FullDescription
        {
            get
            {
                return Accession + "|" + Name + "|" + FullName;
            }
        }

        public string Name { get; private set; }

        public string FullName { get; private set; }

        public bool IsContaminant { get; set; }

        #endregion Public Properties

        #region Public Indexers

        public char this[int zeroBasedIndex]
        {
            get
            {
                return BaseSequence[zeroBasedIndex];
            }
        }

        #endregion Public Indexers

    }
}