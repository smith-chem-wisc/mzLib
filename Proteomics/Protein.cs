using System;
using System.Collections.Generic;

namespace Proteomics
{
    public class Protein
    {
        #region Public Constructors

        public Protein(string sequence, string accession, List<Tuple<string, string>> gene_names = null, IDictionary<int, List<Modification>> oneBasedModifications = null, List<ProteolysisProduct> proteolysisProducts = null, string name = null, string full_name = null, bool isDecoy = false, bool isContaminant = false, List<DatabaseReference> databaseReferences = null, List<SequenceVariation> sequenceVariations = null, List<DisulfideBond> disulfideBonds = null, string databaseFilePath = null)
        {
            // Mandatory
            BaseSequence = sequence;
            Accession = accession;

            Name = name;
            FullName = full_name;
            IsDecoy = isDecoy;
            IsContaminant = isContaminant;
            DatabaseFilePath = databaseFilePath;

            GeneNames = gene_names ?? new List<Tuple<string, string>>();
            ProteolysisProducts = proteolysisProducts ?? new List<ProteolysisProduct>();
            SequenceVariations = sequenceVariations ?? new List<SequenceVariation>();
            OneBasedPossibleLocalizedModifications = oneBasedModifications ?? new Dictionary<int, List<Modification>>();
            DatabaseReferences = databaseReferences ?? new List<DatabaseReference>();
            DisulfideBonds = disulfideBonds ?? new List<DisulfideBond>();
        }

        #endregion Public Constructors

        #region Public Properties

        public IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; }

        /// <summary>
        /// The list of gene names consists of tuples, where Item1 is the type of gene name, and Item2 is the name. There may be many genes and names of a certain type produced when reading an XML protein database.
        /// </summary>
        public IEnumerable<Tuple<string, string>> GeneNames { get; }

        public string Accession { get; }
        public string BaseSequence { get; }
        public bool IsDecoy { get; }
        public IEnumerable<SequenceVariation> SequenceVariations { get; }
        public IEnumerable<DisulfideBond> DisulfideBonds { get; }
        public IEnumerable<ProteolysisProduct> ProteolysisProducts { get; }
        public IEnumerable<DatabaseReference> DatabaseReferences { get; }
        public string DatabaseFilePath { get; }

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

        public string Name { get; }

        public string FullName { get; }

        public bool IsContaminant { get; }

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