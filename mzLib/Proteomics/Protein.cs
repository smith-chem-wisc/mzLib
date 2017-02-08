using System.Collections.Generic;

namespace Proteomics
{
    public class Protein
    {
        #region Private Fields

        private string fullDescription;

        #endregion Private Fields

        #region Public Constructors

        public Protein(string sequence, string accession, Dictionary<int, HashSet<BaseModification>> oneBasedModifications, int[] oneBasedBeginPositions, int[] oneBasedEndPositions, string[] bigPeptideTypes, string name, string full_name, int offset, bool isDecoy, bool isContaminant)
        {
            BaseSequence = sequence;
            Accession = accession;
            OneBasedPossibleLocalizedModifications = oneBasedModifications;
            OneBasedBeginPositions = oneBasedBeginPositions;
            OneBasedEndPositions = oneBasedEndPositions;
            BigPeptideTypes = bigPeptideTypes;
            Name = name;
            FullName = full_name;
            Offset = offset;
            IsDecoy = isDecoy;
            IsContaminant = isContaminant;
        }

        #endregion Public Constructors

        #region Public Properties

        public int[] OneBasedBeginPositions { get; private set; }
        public int[] OneBasedEndPositions { get; private set; }
        public string[] BigPeptideTypes { get; private set; }
        public Dictionary<int, HashSet<BaseModification>> OneBasedPossibleLocalizedModifications { get; private set; }
        public string Accession { get; private set; }
        public string BaseSequence { get; private set; }
        public bool IsDecoy { get; private set; }

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
                if (fullDescription == null)
                {
                    fullDescription = Accession + "|" + Name + "|" + FullName;
                }
                return fullDescription;
            }
        }

        public string Name { get; private set; }

        public string FullName { get; private set; }

        public int Offset { get; private set; }
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