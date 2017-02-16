using System.Collections.Generic;

namespace Proteomics
{
    public class Protein
    {

        #region Public Constructors

        public Protein(string sequence, string accession, IDictionary<int, List<Modification>> oneBasedModifications, int?[] oneBasedBeginPositions, int?[] oneBasedEndPositions, string[] bigPeptideTypes, string name, string full_name, bool isDecoy, bool isContaminant, List<GoTerm> goTerms)
        {
            BaseSequence = sequence;
            Accession = accession;
            OneBasedPossibleLocalizedModifications = oneBasedModifications;
            OneBasedBeginPositions = oneBasedBeginPositions;
            OneBasedEndPositions = oneBasedEndPositions;
            BigPeptideTypes = bigPeptideTypes;
            Name = name;
            FullName = full_name;
            IsDecoy = isDecoy;
            IsContaminant = isContaminant;
            GoTerms = goTerms;
        }

        #endregion Public Constructors

        #region Public Properties

        public int?[] OneBasedBeginPositions { get; private set; }
        public int?[] OneBasedEndPositions { get; private set; }
        public string[] BigPeptideTypes { get; private set; }
        public IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; private set; }
        public string Accession { get; private set; }
        public string BaseSequence { get; private set; }
        public bool IsDecoy { get; private set; }
        public List<GoTerm> GoTerms { get; private set; }

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