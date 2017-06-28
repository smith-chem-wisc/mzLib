using System;
using System.Collections.Generic;

namespace Proteomics
{
    public class Protein
    {

        #region Public Constructors

        public Protein(string sequence, string accession, IEnumerable<Tuple<string, string>> gene_names, IDictionary<int, List<Modification>> oneBasedModifications, int?[] oneBasedBeginPositionsForProteolysisProducts, int?[] oneBasedEndPositionsForProteolysisProducts, string[] oneBasedProteolysisProductsTypes, string name, string full_name, bool isDecoy, bool isContaminant, IEnumerable<DatabaseReference> databaseReferences, IEnumerable<SequenceVariation> sequenceVariations, IEnumerable<DisulfideBond> disulfideBonds)
        {
            BaseSequence = sequence;
            Accession = accession;
            Name = name;
            FullName = full_name;
            IsDecoy = isDecoy;
            IsContaminant = isContaminant;
            GeneNames = gene_names;
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
            SequenceVariations = sequenceVariations;
            DisulfideBonds = disulfideBonds;
            OneBasedPossibleLocalizedModifications = oneBasedModifications;
            DatabaseReferences = databaseReferences;
        }

        #endregion Public Constructors

        #region Public Properties

        public IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; private set; }

        /// <summary>
        /// The list of gene names consists of tuples, where Item1 is the type of gene name, and Item2 is the name. There may be many genes and names of a certain type produced when reading an XML protein database.
        /// </summary>
        public IEnumerable<Tuple<string, string>> GeneNames { get; private set; }

        public string Accession { get; private set; }
        public string BaseSequence { get; private set; }
        public bool IsDecoy { get; private set; }
        public IEnumerable<SequenceVariation> SequenceVariations { get; private set; }
        public IEnumerable<DisulfideBond> DisulfideBonds { get; private set; }
        public IEnumerable<ProteolysisProduct> ProteolysisProducts { get; private set; }
        public IEnumerable<DatabaseReference> DatabaseReferences { get; private set; }

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