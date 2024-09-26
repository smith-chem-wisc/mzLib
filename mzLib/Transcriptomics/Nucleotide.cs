using Chemistry;
using System.Globalization;

namespace Transcriptomics
{
    /// <summary>
    /// Class representing nucleotides
    /// Static members are set of for common bases that are universal and instantiated in the static constructor
    /// Novel bases can be created trough the AddResidue() method
    ///
    /// Future directions could further separate out the sugar and phosphate groups to account for modifications
    /// </summary>
    public class Nucleotide : INucleotide, IEquatable<Nucleotide>
    {
        #region Common Bases

        // RNA
        public static Nucleotide AdenineBase { get; private set; }
        public static Nucleotide CytosineBase { get; private set; }
        public static Nucleotide GuanineBase { get; private set; }
        public static Nucleotide UracilBase { get; private set; }
        public static Nucleotide PseudoUracilBase { get; private set; }

        // DNA
        public static Nucleotide DeoxyAdenineBase { get; private set; }
        public static Nucleotide DeoxyCytosineBase { get; private set; }
        public static Nucleotide DeoxyGuanineBase { get; private set; }
        public static Nucleotide DeoxyThymineBase { get; private set; }

        #endregion

        #region Static Instantiation of Common Known Bases

        /// <summary>
        /// A dictionary of all known residues, can be searched via one letter, three letter, or the name of the residue
        /// </summary>
        internal static readonly Dictionary<string, Nucleotide> AllKnownResidues;

        internal static readonly Nucleotide[] ResiduesByLetter;

        /// <summary>
        /// Constructs all known nucleic acids
        /// </summary>
        static Nucleotide()
        {

            AllKnownResidues = new Dictionary<string, Nucleotide>(66);
            ResiduesByLetter = new Nucleotide['z' + 1]; //Make it big enough for all the Upper and Lower characters

            // actual base chemical formula after bonding with the sugar
            // the sugar and phosphate has a chemical formula of C5H8O6P1
            AdenineBase = AddResidue("Adenine", 'A', "Ade", "C5H4N5");
            CytosineBase = AddResidue("Cytosine", 'C', "Cyt", "C4H4N3O1");
            GuanineBase = AddResidue("Guanine", 'G', "Gua", "C5H4N5O1");
            UracilBase = AddResidue("Uracil", 'U', "Ura", "C4H3N2O2");
            PseudoUracilBase = AddResidue("PseudoUracil", 'Y', "Psu", "C4H3N2O2"); // Y was choosen for pseudouridine due to it commonly being represented by Psi

            // DNA bases which have the same mass as the ones above
            // however, naming to deoxy- to distinguish DNA nucleotide mass calculation from RNA
            DeoxyAdenineBase = AddResidue("DeoxyAdenine", 'B', "dAde", "C5H4N5");
            DeoxyCytosineBase = AddResidue("DeoxyCytosine", 'D', "dCyt", "C4H4N3O1");
            DeoxyGuanineBase = AddResidue("DeoxyGuanine", 'H', "dGua", "C5H4N5O1");
            DeoxyThymineBase = AddResidue("DeoxyThymine", 'V', "dThy", "C5H5N2O2");
        }

        /// <summary>
        /// Adds residue to static internal dictionary of all known residues
        /// </summary>
        /// <param name="residue"></param>
        private static void AddResidueToDictionary(Nucleotide residue)
        {
            AllKnownResidues.Add(residue.Letter.ToString(CultureInfo.InvariantCulture), residue);
            AllKnownResidues.Add(residue.Name, residue);
            AllKnownResidues.Add(residue.Symbol, residue);
            ResiduesByLetter[residue.Letter] = residue;
            ResiduesByLetter[Char.ToLower(residue.Letter)] = residue;
        }

        #endregion

        #region Properties

        /// <summary>
        /// Chemical formula of the sugar less OH 
        /// </summary>
        private ChemicalFormula _sugarAndPhosphate;

        public string Name { get; private set; }
        public double MonoisotopicMass { get; private set; }
        public char Letter { get; private set; }
        public string Symbol { get; private set; }

        /// <summary>
        /// Chemical formula of the base plus the sugar and phosphate less H20
        /// </summary>
        public ChemicalFormula ThisChemicalFormula { get; private set; }

        /// <summary>
        /// Chemical formula of the base less one H
        /// </summary>
        public ChemicalFormula BaseChemicalFormula { get; private set; }

        /// <summary>
        /// Chemical formula of the base plus the sugar
        /// Used to calculate the mass of a modification from the modomics database
        /// </summary>
        public ChemicalFormula NucleosideChemicalFormula { get; private set; }

        #endregion

        #region Constructors

        internal Nucleotide(string name, char oneLetterAbbreviation, string threeLetterAbbreviation, string baseChemicalFormula)
            : this(name, oneLetterAbbreviation, threeLetterAbbreviation, ChemicalFormula.ParseFormula(baseChemicalFormula))
        {
        }

        public Nucleotide(string name, char oneLetterAbbreviation, string threeLetterAbbreviation,
            ChemicalFormula baseChemicalFormula)
        {
            Name = name;
            Letter = oneLetterAbbreviation;
            Symbol = threeLetterAbbreviation;
            ThisChemicalFormula = new ChemicalFormula(baseChemicalFormula);
            BaseChemicalFormula = baseChemicalFormula;

            // calculation for monoisoptic mass of DNA and RNA
            if (Name.Equals("DeoxyAdenine") || Name.Equals("DeoxyCytosine") || Name.Equals("DeoxyGuanine") || Name.Equals("DeoxyThymine"))
            {
                // DNA sugar phosphate backbone (one less oxygen than RNA)
                _sugarAndPhosphate = ChemicalFormula.ParseFormula("C5H8O5P1");
            }
            else
            {
                // RNA sugar phosphate backbone
                _sugarAndPhosphate = ChemicalFormula.ParseFormula("C5H8O6P1");
            }

            ThisChemicalFormula.Add(_sugarAndPhosphate);
            MonoisotopicMass = ThisChemicalFormula.MonoisotopicMass;

            // subtracting one phospho (H03P) and adding one water (H2O) to get the nucleoside formula
            // must add water because base and sugar are already calculated assuming loss of H20 from polymerization reaction
            NucleosideChemicalFormula = ThisChemicalFormula - ChemicalFormula.ParseFormula("H-1O2P1");
        }

        #endregion

        #region Static Methods

        /// <summary>
        /// Adds residue to AllKnownResidues Dictionary
        /// </summary>
        /// <param name="name"></param>
        /// <param name="oneLetterAbbreviation"></param>
        /// <param name="threeLetterAbbreviation"></param>
        /// <param name="chemicalFormula"></param>
        /// <returns></returns>
        public static Nucleotide AddResidue(string name, char oneLetterAbbreviation, string threeLetterAbbreviation, string chemicalFormula)
        {
            var residue = new Nucleotide(name, oneLetterAbbreviation, threeLetterAbbreviation, chemicalFormula);
            AddResidueToDictionary(residue);
            return residue;
        }

        /// <summary>
        /// Get the residue based on the residues' symbol
        /// </summary>
        /// <param name="symbol"></param>
        /// <returns></returns>
        public static Nucleotide GetResidue(string symbol)
        {
            return symbol.Length == 1 ? ResiduesByLetter[symbol[0]] : AllKnownResidues[symbol];
        }

        /// <summary>
        /// Gets the residue based on the residue's one-character symbol
        /// </summary>
        /// <param name="letter"></param>
        /// <returns></returns>
        public static Nucleotide GetResidue(char letter)
        {
            return ResiduesByLetter[letter];
        }

        /// <summary>
        /// Try to get the residue based upon the residue's one-character symbol
        /// </summary>
        /// <param name="letter"></param>
        /// <param name="residue"></param>
        /// <returns></returns>
        public static bool TryGetResidue(char letter, out Nucleotide residue)
        {
            residue = null;
            if (letter > 'z' || letter < 0)
                return false;
            residue = ResiduesByLetter[letter];
            return residue != null;
        }

        /// <summary>
        /// Try to get the residue based upon the residues' symbol
        /// </summary>
        /// <param name="symbol"></param>
        /// <param name="residue"></param>
        /// <returns></returns>
        public static bool TryGetResidue(string symbol, out Nucleotide residue)
        {
            return AllKnownResidues.TryGetValue(symbol, out residue);
        }

        #endregion

        #region Interface Implementations and Overrides

        public override string ToString()
        {
            return $"{Letter} {Symbol} ({Name})";
        }

        public bool Equals(Nucleotide? other)
        {
            if (ReferenceEquals(null, other)) return false;
            if (ReferenceEquals(this, other)) return true;
            return Name == other.Name && ThisChemicalFormula.Equals(other.ThisChemicalFormula)
                                      && Letter == other.Letter
                                      && Symbol == other.Symbol;
        }

        public override bool Equals(object? obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != this.GetType()) return false;
            return Equals((Nucleotide)obj);
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(Name, ThisChemicalFormula, Letter, Symbol);
        }

        #endregion
    }
}
