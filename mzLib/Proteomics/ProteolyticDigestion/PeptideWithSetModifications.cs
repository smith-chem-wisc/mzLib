 using Chemistry;
using MassSpectrometry;
using Proteomics.AminoAcidPolymer;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Proteomics.ProteolyticDigestion
{
    [Serializable]
    public class PeptideWithSetModifications : ProteolyticPeptide
    {
        public string FullSequence { get; private set; } //sequence with modifications
        public readonly int NumFixedMods;
        // Parameter to store a hash code corresponding to a Decoy or a Target peptide
        // If the peptide in question is a decoy, this pairs it to the target it was generated from
        // If the peptide in question is a target, this pairs it to its corresponding decoy 
        public int? PairedTargetDecoyHash { get; private set; }
        /// <summary>
        /// Dictionary of modifications on the peptide. The N terminus is index 1.
        /// The key indicates which residue modification is on (with 1 being N terminus).
        /// </summary>
        [NonSerialized] private Dictionary<int, Modification> _allModsOneIsNterminus; //we currently only allow one mod per position
        [NonSerialized] private bool? _hasChemicalFormulas;
        [NonSerialized] private string _sequenceWithChemicalFormulas;
        [NonSerialized] private double? _monoisotopicMass;
        [NonSerialized] private double? _mostAbundantMonoisotopicMass;
        [NonSerialized] private ChemicalFormula _fullChemicalFormula;
        [NonSerialized] private DigestionParams _digestionParams;
        private static readonly double WaterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private readonly string ProteinAccession; // used to get protein object after deserialization
        /// <summary>
        /// Creates a PeptideWithSetModifications object from a protein. Used when a Protein is digested.
        /// </summary>
        public PeptideWithSetModifications(Protein protein, DigestionParams digestionParams, int oneBasedStartResidueInProtein,
            int oneBasedEndResidueInProtein, CleavageSpecificity cleavageSpecificity, string peptideDescription, int missedCleavages,
           Dictionary<int, Modification> allModsOneIsNterminus, int numFixedMods, string baseSequence = null, int? pairedTargetDecoyHash = null)
           : base(protein, oneBasedStartResidueInProtein, oneBasedEndResidueInProtein, missedCleavages, cleavageSpecificity, peptideDescription, baseSequence)
        {
            _allModsOneIsNterminus = allModsOneIsNterminus;
            NumFixedMods = numFixedMods;
            _digestionParams = digestionParams;
            DetermineFullSequence();
            ProteinAccession = protein.Accession;
            UpdateCleavageSpecificity();
            PairedTargetDecoyHash = pairedTargetDecoyHash; // Added PairedTargetDecoyHash as a nullable integer
        }

        /// <summary>
        /// Creates a PeptideWithSetModifications object from a sequence string.
        /// Useful for reading in MetaMorpheus search engine output into mzLib objects.
        /// </summary>
        public PeptideWithSetModifications(string sequence, Dictionary<string, Modification> allKnownMods, int numFixedMods = 0,
            DigestionParams digestionParams = null, Protein p = null, int oneBasedStartResidueInProtein = int.MinValue,
            int oneBasedEndResidueInProtein = int.MinValue, int missedCleavages = int.MinValue,
            CleavageSpecificity cleavageSpecificity = CleavageSpecificity.Full, string peptideDescription = null, int? pairedTargetDecoyHash = null)
            : base(p, oneBasedStartResidueInProtein, oneBasedEndResidueInProtein, missedCleavages, cleavageSpecificity, peptideDescription)
        {
            if (sequence.Contains("|"))
            {
                throw new MzLibUtil.MzLibException("Ambiguous peptide cannot be parsed from string: " + sequence);
            }

            FullSequence = sequence;
            _baseSequence = GetBaseSequenceFromFullSequence(sequence);
            GetModsAfterDeserialization(allKnownMods);
            NumFixedMods = numFixedMods;
            _digestionParams = digestionParams;
            PairedTargetDecoyHash = pairedTargetDecoyHash; // Added PairedTargetDecoyHash as a nullable integer

            if (p != null)
            {
                ProteinAccession = p.Accession;
            }
        }

        public DigestionParams DigestionParams
        {
            get { return _digestionParams; }
        }

        public Dictionary<int, Modification> AllModsOneIsNterminus
        {
            get { return _allModsOneIsNterminus; }
        }

        public int NumMods
        {
            get { return AllModsOneIsNterminus.Count; }
        }

        public int NumVariableMods
        {
            get { return NumMods - NumFixedMods; }
        }

        public double MonoisotopicMass
        {
            get
            {
                if (!_monoisotopicMass.HasValue)
                {
                    double monoMass = WaterMonoisotopicMass;

                    foreach (var mod in AllModsOneIsNterminus.Values)
                    {
                        monoMass += mod.MonoisotopicMass.Value;
                    }
                    monoMass += BaseSequence.Sum(b => Residue.ResidueMonoisotopicMass[b]);

                    _monoisotopicMass = monoMass;
                }
                return (double)ClassExtensions.RoundedDouble(_monoisotopicMass.Value);
            }

        }
        
        public ChemicalFormula FullChemicalFormula
        {
            get
            {
                ChemicalFormula fullChemicalFormula = new Proteomics.AminoAcidPolymer.Peptide(BaseSequence).GetChemicalFormula();
                foreach (var mod in AllModsOneIsNterminus.Values)
                {
                    if (mod.ChemicalFormula != null)
                    {
                        fullChemicalFormula.Add(mod.ChemicalFormula);
                    }
                    else
                    {
                        fullChemicalFormula = null;
                        break;
                    }
                }
                
                _fullChemicalFormula = fullChemicalFormula;
                return _fullChemicalFormula;
            }
        }

        public double MostAbundantMonoisotopicMass
        {
            get
            {
                if (!_mostAbundantMonoisotopicMass.HasValue)
                {
                    IsotopicDistribution dist = IsotopicDistribution.GetDistribution(this.FullChemicalFormula);
                    double maxIntensity = dist.Intensities.Max();
                    _mostAbundantMonoisotopicMass = (double)ClassExtensions.RoundedDouble(dist.Masses.ToList()[dist.Intensities.ToList().IndexOf(maxIntensity)]);
                }
                return (double)ClassExtensions.RoundedDouble(_mostAbundantMonoisotopicMass.Value);
            }

        }

        public string SequenceWithChemicalFormulas
        {
            get
            {
                if (!_hasChemicalFormulas.HasValue)
                {
                    _hasChemicalFormulas = true;
                    var subsequence = new StringBuilder();

                    // variable modification on peptide N-terminus
                    if (AllModsOneIsNterminus.TryGetValue(1, out Modification pep_n_term_variable_mod))
                    {
                        if (pep_n_term_variable_mod is Modification jj)
                        {
                            subsequence.Append('[' + jj.ChemicalFormula.Formula + ']');
                        }
                        else
                        {
                            return null;
                        }
                    }

                    for (int r = 0; r < Length; r++)
                    {
                        subsequence.Append(this[r]);
                        // variable modification on this residue
                        if (AllModsOneIsNterminus.TryGetValue(r + 2, out Modification residue_variable_mod))
                        {
                            if (residue_variable_mod is Modification jj)
                            {
                                subsequence.Append('[' + jj.ChemicalFormula.Formula + ']');
                            }
                            else
                            {
                                return null;
                            }
                        }
                    }

                    // variable modification on peptide C-terminus
                    if (AllModsOneIsNterminus.TryGetValue(Length + 2, out Modification pep_c_term_variable_mod))
                    {
                        if (pep_c_term_variable_mod is Modification jj)
                        {
                            subsequence.Append('[' + jj.ChemicalFormula.Formula + ']');
                        }
                        else
                        {
                            return null;
                        }
                    }

                    _sequenceWithChemicalFormulas = subsequence.ToString();
                }
                return _sequenceWithChemicalFormulas;
            }
        }

        /// <summary>
        /// Generates theoretical fragments for given dissociation type for this peptide. 
        /// The "products" parameter is filled with these fragments.
        /// </summary>
        public void Fragment(DissociationType dissociationType, FragmentationTerminus fragmentationTerminus, List<Product> products)
        {
            // This code is specifically written to be memory- and CPU -efficient because it is 
            // called millions of times for a typical search (i.e., at least once per peptide). 
            // If you modify this code, BE VERY CAREFUL about allocating new memory, especially 
            // for new collections. This code also deliberately avoids using "yield return", again
            // for performance reasons. Be sure to benchmark any changes with a parallelized 
            // fragmentation of every peptide in a database (i.e., test for speed decreases and 
            // memory issues).

            products.Clear();

            var massCaps = DissociationTypeCollection.GetNAndCTerminalMassShiftsForDissociationType(dissociationType);

            double cTermMass = 0;
            double nTermMass = 0;

            List<ProductType> nTermProductTypes = DissociationTypeCollection.GetTerminusSpecificProductTypesFromDissociation(dissociationType, FragmentationTerminus.N);
            List<ProductType> cTermProductTypes = DissociationTypeCollection.GetTerminusSpecificProductTypesFromDissociation(dissociationType, FragmentationTerminus.C);

            bool calculateNTermFragments = fragmentationTerminus == FragmentationTerminus.N
                || fragmentationTerminus == FragmentationTerminus.Both;

            bool calculateCTermFragments = fragmentationTerminus == FragmentationTerminus.C
                || fragmentationTerminus == FragmentationTerminus.Both;

            //From http://www.matrixscience.com/help/fragmentation_help.html
            //Low Energy CID -- In low energy CID(i.e.collision induced dissociation in a triple quadrupole or an ion trap) a peptide carrying a positive charge fragments mainly along its backbone, 
            //generating predominantly b and y ions. In addition, for fragments containing RKNQ, peaks are seen for ions that have lost ammonia (-17 Da) denoted a*, b* and y*. For fragments containing 
            //STED, loss of water(-18 Da) is denoted a°, b° and y°. Satellite ions from side chain cleavage are not observed.
            bool haveSeenNTermDegreeIon = false;
            bool haveSeenNTermStarIon = false;
            bool haveSeenCTermDegreeIon = false;
            bool haveSeenCTermStarIon = false;

            // these two collections keep track of the neutral losses observed so far on the n-term or c-term.
            // they are apparently necessary, but allocating memory for collections in this function results in
            // inefficient memory usage and thus frequent garbage collection. 
            // TODO: If you can think of a way to remove these collections and still maintain correct 
            // fragmentation, please do so.
            HashSet<double> nTermNeutralLosses = null;
            HashSet<double> cTermNeutralLosses = null;

            // n-terminus mod
            if (calculateNTermFragments)
            {
                if (AllModsOneIsNterminus.TryGetValue(1, out Modification mod))
                {
                    nTermMass += mod.MonoisotopicMass.Value;

                    // n-term mod neutral loss
                    nTermNeutralLosses = AddNeutralLossesFromMods(mod, nTermNeutralLosses, dissociationType);
                }
            }

            // c-terminus mod
            if (calculateCTermFragments)
            {
                if (AllModsOneIsNterminus.TryGetValue(BaseSequence.Length + 2, out Modification mod))
                {
                    cTermMass += mod.MonoisotopicMass.Value;

                    // c-term mod neutral loss
                    cTermNeutralLosses = AddNeutralLossesFromMods(mod, cTermNeutralLosses, dissociationType);
                }
            }

            for (int r = 0; r < BaseSequence.Length - 1; r++)
            {
                // n-term fragments
                if (calculateNTermFragments)
                {
                    char nTermResidue = BaseSequence[r];

                    // get n-term residue mass
                    if (Residue.TryGetResidue(nTermResidue, out Residue residue))
                    {
                        nTermMass += residue.MonoisotopicMass;
                    }
                    else
                    {
                        nTermMass = double.NaN;
                    }

                    // add side-chain mod
                    if (AllModsOneIsNterminus.TryGetValue(r + 2, out Modification mod))
                    {
                        nTermMass += mod.MonoisotopicMass.Value;
                    }

                    // handle star and degree ions for low-res CID
                    if (dissociationType == DissociationType.LowCID)
                    {
                        if (nTermResidue == 'R' || nTermResidue == 'K' || nTermResidue == 'N' || nTermResidue == 'Q')
                        {
                            haveSeenNTermStarIon = true;
                        }

                        if (nTermResidue == 'S' || nTermResidue == 'T' || nTermResidue == 'E' || nTermResidue == 'D')
                        {
                            haveSeenNTermDegreeIon = true;
                        }
                    }

                    // skip first N-terminal fragment (b1, aDegree1, ...) for CID
                    if (r == 0 && (dissociationType == DissociationType.CID || dissociationType == DissociationType.LowCID))
                    {
                        goto CTerminusFragments;
                    }

                    // generate products
                    for (int i = 0; i < nTermProductTypes.Count; i++)
                    {
                        if (dissociationType == DissociationType.LowCID)
                        {
                            if (!haveSeenNTermStarIon && (nTermProductTypes[i] == ProductType.aStar || nTermProductTypes[i] == ProductType.bAmmoniaLoss))
                            {
                                continue;
                            }

                            if (!haveSeenNTermDegreeIon && (nTermProductTypes[i] == ProductType.aDegree || nTermProductTypes[i] == ProductType.bWaterLoss))
                            {
                                continue;
                            }
                        }

                        products.Add(new Product(
                            nTermProductTypes[i],
                            FragmentationTerminus.N,
                            nTermMass + massCaps.Item1[i],
                            r + 1,
                            r + 1,
                            0));

                        nTermNeutralLosses = AddNeutralLossesFromMods(mod, nTermNeutralLosses, dissociationType);

                        if (nTermNeutralLosses != null)
                        {
                            foreach (double neutralLoss in nTermNeutralLosses)
                            {
                                products.Add(new Product(
                                    nTermProductTypes[i],
                                    FragmentationTerminus.N,
                                    nTermMass + massCaps.Item1[i] - neutralLoss,
                                    r + 1,
                                    r + 1,
                                    neutralLoss));
                            }
                        }
                    }
                }

            // c-term fragments
            CTerminusFragments:
                if (calculateCTermFragments)
                {
                    char cTermResidue = BaseSequence[BaseSequence.Length - r - 1];

                    // get c-term residue mass
                    if (Residue.TryGetResidue(cTermResidue, out Residue residue))
                    {
                        cTermMass += residue.MonoisotopicMass;
                    }
                    else
                    {
                        cTermMass = double.NaN;
                    }

                    // add side-chain mod
                    if (AllModsOneIsNterminus.TryGetValue(BaseSequence.Length - r + 1, out Modification mod))
                    {
                        cTermMass += mod.MonoisotopicMass.Value;
                    }

                    // handle star and degree ions for low-res CID
                    if (dissociationType == DissociationType.LowCID)
                    {
                        if (cTermResidue == 'R' || cTermResidue == 'K' || cTermResidue == 'N' || cTermResidue == 'Q')
                        {
                            haveSeenCTermStarIon = true;
                        }

                        if (cTermResidue == 'S' || cTermResidue == 'T' || cTermResidue == 'E' || cTermResidue == 'D')
                        {
                            haveSeenCTermDegreeIon = true;
                        }
                    }

                    // generate products
                    for (int i = 0; i < cTermProductTypes.Count; i++)
                    {
                        // skip zDot ions for proline residues for ETD/ECD/EThcD
                        if (cTermResidue == 'P'
                            && (dissociationType == DissociationType.ECD || dissociationType == DissociationType.ETD || dissociationType == DissociationType.EThcD)
                            && cTermProductTypes[i] == ProductType.zDot)
                        {
                            continue;
                        }

                        if (dissociationType == DissociationType.LowCID)
                        {
                            if (!haveSeenCTermStarIon && cTermProductTypes[i] == ProductType.yAmmoniaLoss)
                            {
                                continue;
                            }

                            if (!haveSeenCTermDegreeIon && cTermProductTypes[i] == ProductType.yWaterLoss)
                            {
                                continue;
                            }
                        }

                        products.Add(new Product(
                            cTermProductTypes[i],
                            FragmentationTerminus.C,
                            cTermMass + massCaps.Item2[i],
                            r + 1,
                            BaseSequence.Length - r,
                            0));

                        cTermNeutralLosses = AddNeutralLossesFromMods(mod, cTermNeutralLosses, dissociationType);

                        if (cTermNeutralLosses != null)
                        {
                            foreach (double neutralLoss in cTermNeutralLosses)
                            {
                                products.Add(new Product(
                                    cTermProductTypes[i],
                                    FragmentationTerminus.C,
                                    cTermMass + massCaps.Item2[i] - neutralLoss,
                                    r + 1,
                                    BaseSequence.Length - r,
                                    neutralLoss));
                            }
                        }
                    }
                }
            }

            // zDot generates one more ion...
            //ETD will cleave between N - C bond.So ETD will remove a NH3 from the N-terminal amino acid, and generate(MH + minus NH3) ion
            if (cTermProductTypes.Contains(ProductType.zDot) && BaseSequence[0] != 'P')
            {
                // get c-term residue mass
                if (Residue.TryGetResidue(BaseSequence[0], out Residue residue))
                {
                    cTermMass += residue.MonoisotopicMass;
                }
                else
                {
                    cTermMass = double.NaN;
                }

                // add side-chain mod
                if (AllModsOneIsNterminus.TryGetValue(2, out Modification mod))
                {
                    cTermMass += mod.MonoisotopicMass.Value;
                }

                // generate zDot product
                products.Add(new Product(
                    ProductType.zDot,
                    FragmentationTerminus.C,
                    cTermMass + DissociationTypeCollection.GetMassShiftFromProductType(ProductType.zDot),
                    BaseSequence.Length,
                    1,
                    0));

                cTermNeutralLosses = AddNeutralLossesFromMods(mod, cTermNeutralLosses, dissociationType);

                if (cTermNeutralLosses != null)
                {
                    foreach (double neutralLoss in cTermNeutralLosses)
                    {
                        products.Add(new Product(
                            ProductType.zDot,
                            FragmentationTerminus.C,
                            cTermMass + DissociationTypeCollection.GetMassShiftFromProductType(ProductType.zDot) - neutralLoss,
                            BaseSequence.Length,
                            1,
                            neutralLoss));
                    }
                }
            }

            foreach (var mod in AllModsOneIsNterminus.Where(p => p.Value.NeutralLosses != null))
            {
                // molecular ion minus neutral losses
                if (mod.Value.NeutralLosses.TryGetValue(dissociationType, out List<double> losses))
                {
                    foreach (double neutralLoss in losses.Where(p => p != 0))
                    {
                        if (neutralLoss != 0)
                        {
                            products.Add(new Product(ProductType.M, FragmentationTerminus.Both, MonoisotopicMass - neutralLoss, 0, 0, neutralLoss));
                        }
                    }
                }

                if (mod.Value.NeutralLosses.TryGetValue(DissociationType.AnyActivationType, out losses))
                {
                    foreach (double neutralLoss in losses.Where(p => p != 0))
                    {
                        if (neutralLoss != 0)
                        {
                            products.Add(new Product(ProductType.M, FragmentationTerminus.Both, MonoisotopicMass - neutralLoss, 0, 0, neutralLoss));
                        }
                    }
                }
            }

            // generate diagnostic ions
            // TODO: this code is memory-efficient but sort of CPU inefficient; it can be further optimized.
            // however, diagnostic ions are fairly rare so it's probably OK for now
            foreach (double diagnosticIon in AllModsOneIsNterminus
                .Where(p => p.Value.DiagnosticIons != null)
                .SelectMany(p => p.Value.DiagnosticIons.Where(v => v.Key == dissociationType || v.Key == DissociationType.AnyActivationType))
                .SelectMany(p => p.Value)
                .Distinct())
            {
                int diagnosticIonLabel = (int)Math.Round(diagnosticIon.ToMz(1), 0);

                // the diagnostic ion is assumed to be annotated in the mod info as the *neutral mass* of the diagnostic ion, not the ionized species
                products.Add(new Product(ProductType.D, FragmentationTerminus.Both, diagnosticIon, diagnosticIonLabel, 0, 0));
            }
        }

        /// <summary>
        /// Generates theoretical internal fragments for given dissociation type for this peptide. 
        /// The "products" parameter is filled with these fragments.
        /// The "minLengthOfFragments" parameter is the minimum number of amino acids for an internal fragment to be included
        /// TODO: Implement neutral losses (e.g. phospho)
        /// TODO: Implement Star/Degree ions from CID
        /// </summary>
        public void FragmentInternally(DissociationType dissociationType, int minLengthOfFragments, List<Product> products)
        {
            products.Clear();

            var massCaps = DissociationTypeCollection.GetNAndCTerminalMassShiftsForDissociationType(dissociationType);

            List<ProductType> nTermProductTypes = DissociationTypeCollection.GetTerminusSpecificProductTypesFromDissociation(dissociationType, FragmentationTerminus.N);
            List<ProductType> cTermProductTypes = DissociationTypeCollection.GetTerminusSpecificProductTypesFromDissociation(dissociationType, FragmentationTerminus.C);

            //foreach start (N-term) index possible
            for (int n = 1; n <= BaseSequence.Length - minLengthOfFragments - 1; n++)
            {
                double fragmentMass = 0;
                //populate with smallest possible fragment (minus 1) from this starting residue
                for (int i = 0; i < minLengthOfFragments - 1; i++)
                {
                    if (Residue.TryGetResidue(BaseSequence[n + i], out Residue residue))
                    {
                        fragmentMass += residue.MonoisotopicMass;

                        // add side-chain mod
                        if (AllModsOneIsNterminus.TryGetValue(n + i + 2, out Modification mod))
                        {
                            fragmentMass += mod.MonoisotopicMass.Value;
                        }
                    }
                    else
                    {
                        fragmentMass = double.NaN;
                    }
                }

                //expand length of fragment, adding each new length as a new fragment ion, until we reach the C1 residue.
                for (int c = n + minLengthOfFragments - 1; c < BaseSequence.Length - 1; c++)
                {
                    if (Residue.TryGetResidue(BaseSequence[c], out Residue residue))
                    {
                        fragmentMass += residue.MonoisotopicMass;
                        // add side-chain mod
                        if (AllModsOneIsNterminus.TryGetValue(c + 2, out Modification mod))
                        {
                            fragmentMass += mod.MonoisotopicMass.Value;
                        }
                        //add new fragment
                        //loop to accomodate EThcD
                        for (int i = 0; i < nTermProductTypes.Count; i++)
                        {
                            double massCap = massCaps.Item1[i];
                            for (int j = 0; j < cTermProductTypes.Count; j++)
                            {
                                double massCap2 = massCaps.Item2[j];
                                //do c, then n terminal ions
                                products.Add(new Product(cTermProductTypes[j], FragmentationTerminus.None, fragmentMass + massCap + massCap2 - WaterMonoisotopicMass,
                                    n + 1, c - n + 1, 0, nTermProductTypes[i], c + 1));
                            }
                        }
                    }
                    else
                    {
                        fragmentMass = double.NaN;
                    }
                }
            }
        }

        public virtual string EssentialSequence(IReadOnlyDictionary<string, int> modstoWritePruned)
        {
            string essentialSequence = BaseSequence;
            if (modstoWritePruned != null)
            {
                var sbsequence = new StringBuilder();

                // variable modification on peptide N-terminus
                if (AllModsOneIsNterminus.TryGetValue(1, out Modification pep_n_term_variable_mod))
                {
                    if (modstoWritePruned.ContainsKey(pep_n_term_variable_mod.ModificationType))
                    {
                        sbsequence.Append('[' + pep_n_term_variable_mod.ModificationType + ":" + pep_n_term_variable_mod.IdWithMotif + ']');
                    }
                }
                for (int r = 0; r < Length; r++)
                {
                    sbsequence.Append(this[r]);
                    // variable modification on this residue
                    if (AllModsOneIsNterminus.TryGetValue(r + 2, out Modification residue_variable_mod))
                    {
                        if (modstoWritePruned.ContainsKey(residue_variable_mod.ModificationType))
                        {
                            sbsequence.Append('[' + residue_variable_mod.ModificationType + ":" + residue_variable_mod.IdWithMotif + ']');
                        }
                    }
                }

                // variable modification on peptide C-terminus
                if (AllModsOneIsNterminus.TryGetValue(Length + 2, out Modification pep_c_term_variable_mod))
                {
                    if (modstoWritePruned.ContainsKey(pep_c_term_variable_mod.ModificationType))
                    {
                        sbsequence.Append('[' + pep_c_term_variable_mod.ModificationType + ":" + pep_c_term_variable_mod.IdWithMotif + ']');
                    }
                }

                essentialSequence = sbsequence.ToString();
            }
            return essentialSequence;
        }

        public PeptideWithSetModifications Localize(int j, double massToLocalize)
        {
            var dictWithLocalizedMass = new Dictionary<int, Modification>(AllModsOneIsNterminus);
            double massOfExistingMod = 0;
            if (dictWithLocalizedMass.TryGetValue(j + 2, out Modification modToReplace))
            {
                massOfExistingMod = (double)modToReplace.MonoisotopicMass;
                dictWithLocalizedMass.Remove(j + 2);
            }

            dictWithLocalizedMass.Add(j + 2, new Modification(_locationRestriction: "Anywhere.", _monoisotopicMass: massToLocalize + massOfExistingMod));

            var peptideWithLocalizedMass = new PeptideWithSetModifications(Protein, _digestionParams, OneBasedStartResidueInProtein, OneBasedEndResidueInProtein,
                CleavageSpecificityForFdrCategory, PeptideDescription, MissedCleavages, dictWithLocalizedMass, NumFixedMods);

            return peptideWithLocalizedMass;
        }

        /// <summary>
        /// Determines whether a peptide includes a splice site
        /// </summary>
        /// <param name="pep"></param>
        /// <param name="site"></param>
        /// <returns></returns>
        public bool IncludesSpliceSite(SpliceSite site)
        {
            return OneBasedStartResidueInProtein <= site.OneBasedBeginPosition && OneBasedEndResidueInProtein >= site.OneBasedEndPosition;
        }

        /// <summary>
        /// Checks if sequence variant and peptide intersect, also checks if the seuqence variatn can be identified whether they intersect
        /// or not (ie if the variant causes a cleavage site generating the peptide). Returns a tuple with item 1 being a bool value
        /// representing if the varaint intersects the peptide and item 2 beign abool that represents if the variatn is identified.
        /// </summary>
        /// <param name="pep"></param>
        /// <param name="appliedVariation"></param>
        /// <returns></returns>
        public (bool intersects, bool identifies) IntersectsAndIdentifiesVariation(SequenceVariation appliedVariation)
        {
            // does it intersect?
            //possible locations for variant start site
            bool VariantStartsBeforePeptide = appliedVariation.OneBasedBeginPosition < OneBasedStartResidueInProtein;
            bool VariantStartsAtPeptideStart = appliedVariation.OneBasedBeginPosition == OneBasedStartResidueInProtein;
            bool VariantStartsInsidePeptide = appliedVariation.OneBasedBeginPosition >= OneBasedStartResidueInProtein && appliedVariation.OneBasedBeginPosition < OneBasedEndResidueInProtein;
            bool VariantStartsAtPeptideEnd = appliedVariation.OneBasedBeginPosition == OneBasedEndResidueInProtein;
            //possibe locations for variant end stite
            bool VariantEndsAtPeptideStart = appliedVariation.OneBasedEndPosition == OneBasedStartResidueInProtein;
            bool VariantEndsInsidePeptide = appliedVariation.OneBasedEndPosition > OneBasedStartResidueInProtein && appliedVariation.OneBasedEndPosition <= OneBasedEndResidueInProtein;
            bool VariantEndsAtPeptideEnd = appliedVariation.OneBasedEndPosition == OneBasedEndResidueInProtein;
            bool VariantEndsAfterPeptide = appliedVariation.OneBasedEndPosition > OneBasedEndResidueInProtein;

            bool intersects = false;
            bool identifies = false;
            //start and end  combinations that lead to variants being intersected by the peptide sequnce
            if (VariantStartsBeforePeptide || VariantStartsAtPeptideStart)
            {
                if (VariantEndsAtPeptideStart || VariantEndsInsidePeptide || VariantEndsAtPeptideEnd || VariantEndsAfterPeptide)
                {
                    intersects = true;
                }
            }
            else if (VariantStartsInsidePeptide)
            {
                if (VariantEndsInsidePeptide || VariantEndsAfterPeptide || VariantEndsAtPeptideEnd)
                {
                    intersects = true;
                }
            }
            else if (VariantStartsAtPeptideEnd)
            {
                if (VariantEndsAfterPeptide || VariantEndsAtPeptideEnd)
                {
                    intersects = true;
                }
            }

            if (intersects == true)
            {
                int lengthDiff = appliedVariation.VariantSequence.Length - appliedVariation.OriginalSequence.Length;
                int intersectOneBasedStart = Math.Max(OneBasedStartResidueInProtein, appliedVariation.OneBasedBeginPosition);
                int intersectOneBasedEnd = Math.Min(OneBasedEndResidueInProtein, appliedVariation.OneBasedEndPosition + lengthDiff);
                int intersectSize = intersectOneBasedEnd - intersectOneBasedStart + 1;

                // if the original sequence within the peptide is shorter or longer than the variant sequence within the peptide, there is a sequence change
                int variantZeroBasedStartInPeptide = intersectOneBasedStart - appliedVariation.OneBasedBeginPosition;
                bool origSeqIsShort = appliedVariation.OriginalSequence.Length - variantZeroBasedStartInPeptide < intersectSize;
                bool origSeqIsLong = appliedVariation.OriginalSequence.Length > intersectSize && OneBasedEndResidueInProtein > intersectOneBasedEnd;
                if (origSeqIsShort || origSeqIsLong)
                {
                    identifies = true;
                }
                else
                {
                    // crosses the entire variant sequence (needed to identify truncations and certain deletions, like KAAAAAAAAA -> K, but also catches synonymous variations A -> A)
                    bool crossesEntireVariant = intersectSize == appliedVariation.VariantSequence.Length;

                    if (crossesEntireVariant == true)
                    {
                        // is the variant sequence intersecting the peptide different than the original sequence?
                        string originalAtIntersect = appliedVariation.OriginalSequence.Substring(intersectOneBasedStart - appliedVariation.OneBasedBeginPosition, intersectSize);
                        string variantAtIntersect = appliedVariation.VariantSequence.Substring(intersectOneBasedStart - appliedVariation.OneBasedBeginPosition, intersectSize);
                        identifies = originalAtIntersect != variantAtIntersect;
                    }
                }
            }
            //checks to see if the variant causes a cleavage event creating the peptide. This is how a variant can be identified without intersecting
            //with the peptide itself
            else
            {
                //We need to account for any variants that occur in the protien prior to the variant in question.
                //This information is used to calculate a scaling factor to calculate the AA that proceeds the peptide seqeunce in the original (variant free) protein
                List<SequenceVariation> VariantsThatAffectPreviousAAPosition = Protein.AppliedSequenceVariations.Where(v => v.OneBasedEndPosition <= OneBasedStartResidueInProtein).ToList();
                int totalLengthDifference = 0;
                foreach (var variant in VariantsThatAffectPreviousAAPosition)
                {
                    totalLengthDifference += variant.VariantSequence.Length - variant.OriginalSequence.Length;
                }

                //need to determine what the cleavage sites are for the protease used (will allow us to determine if new cleavage sites were made by variant)
                List<DigestionMotif> proteasesCleavageSites = DigestionParams.Protease.DigestionMotifs;
                //if the variant ends the AA before the peptide starts then it may have caused c-terminal cleavage
                //see if the protease used for digestion has C-terminal cleavage sites
                List<string> cTerminalResidue = proteasesCleavageSites.Where(dm => dm.CutIndex == 1).Select(d => d.InducingCleavage).ToList();

                if (appliedVariation.OneBasedEndPosition == (OneBasedStartResidueInProtein - 1))
                {
                    if (cTerminalResidue.Count > 0)
                    {
                        // get the AA that proceeds the peptide from the variant protein (AKA the last AA in the variant)
                        PeptideWithSetModifications previousAA_Variant = new PeptideWithSetModifications(Protein, DigestionParams, OneBasedStartResidueInProtein - 1, OneBasedStartResidueInProtein - 1, CleavageSpecificity.Full, "full", 0, AllModsOneIsNterminus, NumFixedMods);

                        // get the AA that proceeds the peptide sequence in the original protein (wihtout any applied variants)
                        PeptideWithSetModifications previousAA_Original = new PeptideWithSetModifications(Protein.NonVariantProtein, DigestionParams, (OneBasedStartResidueInProtein - 1) - totalLengthDifference, (OneBasedStartResidueInProtein - 1) - totalLengthDifference, CleavageSpecificity.Full, "full", 0, AllModsOneIsNterminus, NumFixedMods);
                        bool newSite = cTerminalResidue.Contains(previousAA_Variant.BaseSequence);
                        bool oldSite = cTerminalResidue.Contains(previousAA_Original.BaseSequence);
                        // if the new AA causes a cleavage event, and that cleavage event would not have occurred without the variant then it is identified
                        if (newSite == true && oldSite == false)
                        {
                            identifies = true;
                        }
                    }
                }
                //if the variant begins the AA after the peptide ends then it may have caused n-terminal cleavage
                else if (appliedVariation.OneBasedBeginPosition == (OneBasedEndResidueInProtein + 1))
                {
                    //see if the protease used for digestion has N-terminal cleavage sites
                    List<string> nTerminalResidue = proteasesCleavageSites.Where(dm => dm.CutIndex == 0).Select(d => d.InducingCleavage).ToList();
                    // stop gain variation can create a peptide this checks for this with cTerminal cleavage proteases
                    if (cTerminalResidue.Count > 0)
                    {
                        if (appliedVariation.VariantSequence == "*")
                        {
                            PeptideWithSetModifications lastAAofPeptide = new PeptideWithSetModifications(Protein, DigestionParams, OneBasedEndResidueInProtein, OneBasedEndResidueInProtein, CleavageSpecificity.Full, "full", 0, AllModsOneIsNterminus, NumFixedMods);
                            bool oldSite = cTerminalResidue.Contains(lastAAofPeptide.BaseSequence);
                            if (oldSite == false)
                            {
                                identifies = true;
                            }
                        }
                    }

                    if (nTerminalResidue.Count > 0)
                    {
                        if (Protein.Length >= OneBasedEndResidueInProtein + 1)
                        {
                            //get the AA that follows the peptide sequence fromt he variant protein (AKA the first AA of the varaint)
                            PeptideWithSetModifications nextAA_Variant = new PeptideWithSetModifications(Protein, DigestionParams, OneBasedEndResidueInProtein + 1, OneBasedEndResidueInProtein + 1, CleavageSpecificity.Full, "full", 0, AllModsOneIsNterminus, NumFixedMods);

                            // checks to make sure the original protein has an amino acid following the peptide (an issue with stop loss variants or variatns that add AA after the previous stop residue)
                            // no else statement because if the peptide end residue was the previous protein stop site, there is no way to truly identify the variant. 
                            // if the peptide were to extend into the stop loss region then the peptide would intesect the variant and this code block would not be triggered.
                            if (Protein.NonVariantProtein.Length >= OneBasedEndResidueInProtein + 1)
                            {
                                // get the AA that follows the peptide sequence in the original protein (without any applied variants)
                                PeptideWithSetModifications nextAA_Original = new PeptideWithSetModifications(Protein.NonVariantProtein, DigestionParams, (OneBasedEndResidueInProtein + 1) - totalLengthDifference, (OneBasedEndResidueInProtein + 1) - totalLengthDifference, CleavageSpecificity.Full, "full", 0, AllModsOneIsNterminus, NumFixedMods);
                                bool newSite = nTerminalResidue.Contains(nextAA_Variant.BaseSequence);
                                bool oldSite = nTerminalResidue.Contains(nextAA_Original.BaseSequence);
                                // if the new AA causes a cleavage event, and that cleavage event would not have occurred without the variant then it is identified
                                if (newSite == true && oldSite == false)
                                {
                                    identifies = true;
                                }
                            }

                        }
                        //for stop gain varations that cause peptide
                        else
                        {
                            // get the AA that follows the peptide sequence in the original protein (without any applied variants)
                            PeptideWithSetModifications nextAA_Original = new PeptideWithSetModifications(Protein.NonVariantProtein, DigestionParams, (OneBasedEndResidueInProtein + 1) - totalLengthDifference, (OneBasedEndResidueInProtein + 1) - totalLengthDifference, CleavageSpecificity.Full, "full", 0, AllModsOneIsNterminus, NumFixedMods);
                            bool oldSite = nTerminalResidue.Contains(nextAA_Original.BaseSequence);
                            // if the new AA causes a cleavage event, and that cleavage event would not have occurred without the variant then it is identified
                            if (oldSite == false)
                            {
                                identifies = true;
                            }
                        }
                    }
                }
            }

            return (intersects, identifies);
        }

        /// <summary>
        /// Makes the string representing a detected sequence variation, including any modifications on a variant amino acid.
        /// takes in the variant as well as the bool value of wheter the peptid eintersects the variant. (this allows for identified
        /// variants that cause the cleavage site for the peptide.
        /// </summary>
        /// <param name="p"></param>
        /// <param name="d"></param>
        /// <returns></returns>
        public string SequenceVariantString(SequenceVariation applied, bool intersects)
        {
            if (intersects == true)
            {
                bool startAtNTerm = applied.OneBasedBeginPosition == 1 && OneBasedStartResidueInProtein == 1;
                bool onlyPeptideStartAtNTerm = OneBasedStartResidueInProtein == 1 && applied.OneBasedBeginPosition != 1;
                int modResidueScale = 0;
                if (startAtNTerm)
                {
                    modResidueScale = 1;
                }
                else if (onlyPeptideStartAtNTerm)
                {
                    modResidueScale = 2;
                }
                else
                {
                    modResidueScale = 3;
                }
                int lengthDiff = applied.VariantSequence.Length - applied.OriginalSequence.Length;
                var modsOnVariantOneIsNTerm = AllModsOneIsNterminus
                    .Where(kv => kv.Key == 1 && applied.OneBasedBeginPosition == 1 || applied.OneBasedBeginPosition <= kv.Key - 2 + OneBasedStartResidueInProtein && kv.Key - 2 + OneBasedStartResidueInProtein <= applied.OneBasedEndPosition)
                    .ToDictionary(kv => kv.Key - applied.OneBasedBeginPosition + (modResidueScale), kv => kv.Value);
                PeptideWithSetModifications variantWithAnyMods = new PeptideWithSetModifications(Protein, DigestionParams, applied.OneBasedBeginPosition == 1 ? applied.OneBasedBeginPosition : applied.OneBasedBeginPosition - 1, applied.OneBasedEndPosition, CleavageSpecificityForFdrCategory, PeptideDescription, MissedCleavages, modsOnVariantOneIsNTerm, NumFixedMods);
                return $"{applied.OriginalSequence}{applied.OneBasedBeginPosition}{variantWithAnyMods.FullSequence.Substring(applied.OneBasedBeginPosition == 1 ? 0 : 1)}";
            }
            //if the variant caused a cleavage site leading the the peptide sequence (variant does not intersect but is identified)
            else
            {
                return $"{applied.OriginalSequence}{ applied.OneBasedBeginPosition}{applied.VariantSequence}";
            }
        }

        /// <summary>
        /// Takes an individual peptideWithSetModifications and determines if applied variations from the protein are found within its length
        /// </summary>
        /// <returns></returns>
        public bool IsVariantPeptide()
        {
            bool identifiedVariant = false;
            if (this.Protein.AppliedSequenceVariations.Count() > 0)
            {
                foreach (var variant in this.Protein.AppliedSequenceVariations)
                {
                    if (this.IntersectsAndIdentifiesVariation(variant).identifies)
                    {
                        identifiedVariant = true;
                        break;
                    }
                }
            }
            return identifiedVariant;
        }

        public override string ToString()
        {
            return FullSequence + string.Join("\t", AllModsOneIsNterminus.Select(m => m.ToString()));
        }

        public override bool Equals(object obj)
        {
            var q = obj as PeptideWithSetModifications;

            if (Protein == null && q.Protein == null)
            {
                return q.FullSequence.Equals(this.FullSequence);
            }

            return q != null
                && q.FullSequence.Equals(this.FullSequence)
                && q.OneBasedStartResidueInProtein == this.OneBasedStartResidueInProtein
                && (q.Protein.Accession == null && this.Protein.Accession == null || q.Protein.Accession.Equals(this.Protein.Accession))
                && q.DigestionParams.Protease.Equals(this.DigestionParams.Protease);
        }

        public override int GetHashCode()
        {
            if (DigestionParams == null)
            {
                return FullSequence.GetHashCode();
            }
            else
            {
                return FullSequence.GetHashCode() + DigestionParams.Protease.GetHashCode();
            }
        }

        /// <summary>
        /// This should be run after deserialization of a PeptideWithSetModifications, in order to set its Protein and Modification objects, which were not serialized
        /// </summary>
        public void SetNonSerializedPeptideInfo(Dictionary<string, Modification> idToMod, Dictionary<string, Protein> accessionToProtein, DigestionParams dp)
        {
            GetModsAfterDeserialization(idToMod);
            GetProteinAfterDeserialization(accessionToProtein);
            _digestionParams = dp;
        }

        private void GetModsAfterDeserialization(Dictionary<string, Modification> idToMod)
        {
            _allModsOneIsNterminus = new Dictionary<int, Modification>();
            int currentModStart = 0;
            int currentModificationLocation = 1;
            bool currentlyReadingMod = false;
            int bracketCount = 0;

            for (int r = 0; r < FullSequence.Length; r++)
            {
                char c = FullSequence[r];
                if (c == '[')
                {
                    currentlyReadingMod = true;
                    if (bracketCount == 0)
                    {
                        currentModStart = r + 1;
                    }
                    bracketCount++;
                }
                else if (c == ']')
                {
                    string modId = null;
                    bracketCount--;
                    if (bracketCount == 0)
                    {
                        try
                        {
                            //remove the beginning section (e.g. "Fixed", "Variable", "Uniprot")
                            string modString = FullSequence.Substring(currentModStart, r - currentModStart);
                            int splitIndex = modString.IndexOf(':');
                            string modType = modString.Substring(0, splitIndex);
                            modId = modString.Substring(splitIndex + 1, modString.Length - splitIndex - 1);
                        }
                        catch (Exception e)
                        {
                            throw new MzLibUtil.MzLibException(
                                "Error while trying to parse string into peptide: " + e.Message);
                        }
                        if (!idToMod.TryGetValue(modId, out Modification mod))
                        {
                            throw new MzLibUtil.MzLibException(
                                "Could not find modification while reading string: " + FullSequence);
                        }
                        if (mod.LocationRestriction.Contains("C-terminal.") && r == FullSequence.Length - 1)
                        {
                            currentModificationLocation = BaseSequence.Length + 2;
                        }
                        _allModsOneIsNterminus.Add(currentModificationLocation, mod);
                        currentlyReadingMod = false;
                    }
                }
                else if (!currentlyReadingMod)
                {
                    currentModificationLocation++;
                }
                //else do nothing
            }
        }

        private void GetProteinAfterDeserialization(Dictionary<string, Protein> idToProtein)
        {
            Protein protein = null;

            if (ProteinAccession != null && !idToProtein.TryGetValue(ProteinAccession, out protein))
            {
                throw new MzLibUtil.MzLibException("Could not find protein accession after deserialization! " + ProteinAccession);
            }

            Protein = protein;
        }

        public static string GetBaseSequenceFromFullSequence(string fullSequence)
        {
            StringBuilder sb = new StringBuilder();
            int bracketCount = 0;
            foreach (char c in fullSequence)
            {
                if (c == '[')
                {
                    bracketCount++;
                }
                else if (c == ']')
                {
                    bracketCount--;
                }
                else if (bracketCount == 0)
                {
                    sb.Append(c);
                }
            }
            return sb.ToString();
        }

        private void DetermineFullSequence()
        {
            var subsequence = new StringBuilder();

            // modification on peptide N-terminus
            if (AllModsOneIsNterminus.TryGetValue(1, out Modification mod))
            {
                subsequence.Append('[' + mod.ModificationType + ":" + mod.IdWithMotif + ']');
            }

            for (int r = 0; r < Length; r++)
            {
                subsequence.Append(this[r]);

                // modification on this residue
                if (AllModsOneIsNterminus.TryGetValue(r + 2, out mod))
                {
                    subsequence.Append('[' + mod.ModificationType + ":" + mod.IdWithMotif + ']');
                }
            }

            // modification on peptide C-terminus
            if (AllModsOneIsNterminus.TryGetValue(Length + 2, out mod))
            {
                subsequence.Append('[' + mod.ModificationType + ":" + mod.IdWithMotif + ']');
            }

            FullSequence = subsequence.ToString();
        }

        private void UpdateCleavageSpecificity()
        {
            if (CleavageSpecificityForFdrCategory == CleavageSpecificity.Unknown)
            {
                CleavageSpecificityForFdrCategory = DigestionParams.SpecificProtease.GetCleavageSpecificity(Protein, OneBasedStartResidueInProtein, OneBasedEndResidueInProtein, DigestionParams.InitiatorMethionineBehavior == InitiatorMethionineBehavior.Retain);
                PeptideDescription = CleavageSpecificityForFdrCategory.ToString();
            }
        }
        
        private HashSet<double> AddNeutralLossesFromMods(Modification mod, HashSet<double> allNeutralLossesSoFar, DissociationType dissociationType)
        {
            // add neutral losses specific to this dissociation type
            if (mod != null
                && mod.NeutralLosses != null
                && mod.NeutralLosses.TryGetValue(dissociationType, out List<double> neutralLossesFromMod))
            {
                foreach (double neutralLoss in neutralLossesFromMod.Where(p => p != 0))
                {
                    if (allNeutralLossesSoFar == null)
                    {
                        allNeutralLossesSoFar = new HashSet<double>();
                    }

                    allNeutralLossesSoFar.Add(neutralLoss);
                }
            }

            // add neutral losses that are generic to any dissociation type
            if (mod != null
                && mod.NeutralLosses != null
                && mod.NeutralLosses.TryGetValue(DissociationType.AnyActivationType, out neutralLossesFromMod))
            {
                foreach (double neutralLoss in neutralLossesFromMod.Where(p => p != 0))
                {
                    if (allNeutralLossesSoFar == null)
                    {
                        allNeutralLossesSoFar = new HashSet<double>();
                    }

                    allNeutralLossesSoFar.Add(neutralLoss);
                }
            }

            return allNeutralLossesSoFar;
        }

        //This function maintains the amino acids associated with the protease motif and reverses all other amino acids.
        //N-terminal modificatons are preserved. Other modifications travel with their respective amino acids. this results
        //in a decoy peptide composed the same amino acids and modifications as the original. 
        //Occasionally, this process results in peptide with exactly the same sequence. Therefore, there is a stop-gap measure
        //the returns the mirror image of the original. N-terminal mods are preserved, but other mods are also reversed. 
        //this should yield a unique decoy for each target sequence.
        //This function also adds a hash code to both the original PeptideWithSetModifications and the decoy 
        //generated by this function pairing the two together by eachother's FullSequence.
        //The original taget peptide is given a hash code corresponding to the decoy's full sequence,
        //and the decoy is given a hash code corresponding to the original target peptide's sequence.
        //This hash code is stored in the PairedTargetDecoyHash parameter of PeptideWithSetModifications.
        public PeptideWithSetModifications GetReverseDecoyFromTarget(int[] revisedAminoAcidOrder)
        {
            Dictionary<int, Modification> newModificationsDictionary = new Dictionary<int, Modification>();
            //Copy N-terminal modifications from target dictionary to decoy dictionary.
            if (this.AllModsOneIsNterminus.ContainsKey(1))
            {
                newModificationsDictionary.Add(1, this.AllModsOneIsNterminus[1]);
            }
            char[] newBase = new char[this.BaseSequence.Length];
            Array.Fill(newBase, '0');
            char[] evaporatingBase = this.BaseSequence.ToCharArray();
            List<DigestionMotif> motifs = this.DigestionParams.Protease.DigestionMotifs;
            if (motifs != null && motifs.Count > 0)
            {
                foreach (var motif in motifs.Where(m => m.InducingCleavage != ""))//check the empty "" for topdown
                {
                    string cleavingMotif = motif.InducingCleavage;
                    List<int> cleavageMotifLocations = new List<int>();

                    for (int i = 0; i < BaseSequence.Length; i++)
                    {
                        bool fits;
                        bool prevents;
                        (fits, prevents) = motif.Fits(BaseSequence, i);

                        if (fits && !prevents)
                        {
                            cleavageMotifLocations.Add(i);
                        }
                    }

                    foreach (int location in cleavageMotifLocations)
                    {
                        char[] motifArray = BaseSequence.Substring(location, cleavingMotif.Length).ToCharArray();

                        for (int i = 0; i < cleavingMotif.Length; i++)
                        {
                            newBase[location + i] = motifArray[i];
                            revisedAminoAcidOrder[location + i] = location + i;//
                            //directly copy mods that were on amino acids in the motif. Those amino acids don't change position.
                            if (this.AllModsOneIsNterminus.ContainsKey(location + i + 2))
                            {
                                newModificationsDictionary.Add(location + i + 2, this.AllModsOneIsNterminus[location + i + 2]);
                            }

                            evaporatingBase[location + i] = '0';//can null a char so i use a number which doesnt' appear in peptide string
                        }
                    }
                }
            }

            // We've kept amino acids in the digestion motif in the same position in the decoy peptide.
            // Now we will fill the remaining open positions in the decoy with the reverse of amino acids from the target.
            // Part to change to scramble
            int fillPosition = 0;
            int extractPosition = this.BaseSequence.Length - 1;
            while (fillPosition < this.BaseSequence.Length && extractPosition >= 0)
            {
                if (evaporatingBase[extractPosition] != '0')
                {
                    while (newBase[fillPosition] != '0')
                    {
                        fillPosition++;
                    }
                    newBase[fillPosition] = evaporatingBase[extractPosition];
                    revisedAminoAcidOrder[fillPosition] = extractPosition;
                    if (this.AllModsOneIsNterminus.ContainsKey(extractPosition + 2))
                    {
                        newModificationsDictionary.Add(fillPosition + 2, this.AllModsOneIsNterminus[extractPosition + 2]);
                    }
                    fillPosition++;
                }
                extractPosition--;
            }

            string newBaseString = new string(newBase);

            var proteinSequence = this.Protein.BaseSequence;
            var aStringBuilder = new StringBuilder(proteinSequence);
            aStringBuilder.Remove(this.OneBasedStartResidueInProtein - 1, this.BaseSequence.Length);
            aStringBuilder.Insert(this.OneBasedStartResidueInProtein - 1, newBaseString);
            proteinSequence = aStringBuilder.ToString();

            Protein decoyProtein = new Protein(proteinSequence, "DECOY_" + this.Protein.Accession, null, new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), null, null, null, true);
            DigestionParams d = this.DigestionParams;

            // Creates a hash code corresponding to the target's sequence
            int targetHash = GetHashCode();
            PeptideWithSetModifications decoyPeptide;
            //Make the "peptideDescription" store the corresponding target's sequence
            if (newBaseString != this.BaseSequence)
            {
                decoyPeptide = new PeptideWithSetModifications(decoyProtein, d, this.OneBasedStartResidueInProtein, this.OneBasedEndResidueInProtein, this.CleavageSpecificityForFdrCategory, this.FullSequence, this.MissedCleavages, newModificationsDictionary, this.NumFixedMods, newBaseString);
                // Sets PairedTargetDecoyHash of the original target peptie to the hash hode of the decoy sequence           
                PairedTargetDecoyHash = decoyPeptide.GetHashCode();
                // Sets PairedTargetDecoyHash of the decoy peptide to the hash code of the target sequence
                decoyPeptide.PairedTargetDecoyHash = targetHash;
                return decoyPeptide;

            }
            else
            {
                //The reverse decoy procedure failed to create a PeptideWithSetModificatons with a different sequence. Therefore,
                //we retrun the mirror image peptide.
                decoyPeptide = this.GetPeptideMirror(revisedAminoAcidOrder);
                PairedTargetDecoyHash = decoyPeptide.GetHashCode();
                decoyPeptide.PairedTargetDecoyHash = targetHash;
                return decoyPeptide;
            }

        }
        /// <summary>
        /// This function generates a decoy peptide from a target by scrambling the target peptide's amino acid sequence
        /// This preserves any digestion motifs and keeps modifications with their amino acids
        /// To help generate only high quality decoys, a homology cutoff of 30 % sequence similarity is used
        /// If after 10 attempts no sufficient decoy is generated, the mirror sequence is returned
        /// </summary>
        /// <param name="revisedAminoAcidOrder">Array to store the new amino acid order in</param>
        /// <param name="maximumHomology">Parameter specifying the homology cutoff to be used</param>
        /// <returns></returns>
        public PeptideWithSetModifications GetScrambledDecoyFromTarget(int[] revisedAminoAcidOrder, double maximumHomology = 0.3)
        {
            Dictionary<int, Modification> newModificationsDictionary = new Dictionary<int, Modification>();
            //Copy N-terminal modifications from target dictionary to decoy dictionary.
            if (this.AllModsOneIsNterminus.ContainsKey(1))
            {
                newModificationsDictionary.Add(1, this.AllModsOneIsNterminus[1]);
            }
            char[] newBase = new char[this.BaseSequence.Length];
            Array.Fill(newBase, '0');
            char[] evaporatingBase = this.BaseSequence.ToCharArray();
            List<DigestionMotif> motifs = this.DigestionParams.Protease.DigestionMotifs;
            if (motifs != null && motifs.Count > 0)
            {
                foreach (var motif in motifs.Where(m => m.InducingCleavage != ""))//check the empty "" for topdown
                {
                    string cleavingMotif = motif.InducingCleavage;
                    List<int> cleavageMotifLocations = new List<int>();

                    for (int i = 0; i < BaseSequence.Length; i++)
                    {
                        bool fits;
                        bool prevents;
                        (fits, prevents) = motif.Fits(BaseSequence, i);

                        if (fits && !prevents)
                        {
                            cleavageMotifLocations.Add(i);
                        }
                    }

                    foreach (int location in cleavageMotifLocations)
                    {
                        char[] motifArray = BaseSequence.Substring(location, cleavingMotif.Length).ToCharArray();
                        
                        for (int i = 0; i < cleavingMotif.Length; i++)
                        {
                            newBase[location + i] = motifArray[i];
                            revisedAminoAcidOrder[location + i] = location + i;
                            //directly copy mods that were on amino acids in the motif. Those amino acids don't change position.
                            if (this.AllModsOneIsNterminus.ContainsKey(location + i + 2))
                            {
                                newModificationsDictionary.Add(location + i + 2, this.AllModsOneIsNterminus[location + i + 2]);
                            }

                            evaporatingBase[location + i] = '0';//can null a char so i use a number which doesnt' appear in peptide string
                        }
                    }
                }
            }

            //We've kept amino acids in the digestion motif in the same position in the decoy peptide.
            //Now we will fill the remaining open positions in the decoy with the scrambled amino acids from the target.
            int extractPosition;
            int fillPosition;
            int residueNumsIndex;
            // Specify seed to ensure that the same decoy sequence is always generated from the target
            Random rand = new(56);
            double percentIdentity = 1;
            int scrambleAttempt = 0;
            int maxScrambles = 10;
            double maxIdentity = maximumHomology;
            int characterCounter;

            while(scrambleAttempt < maxScrambles && percentIdentity > maxIdentity)
            {
                // Copies the newModificationsDictionary for the scramble attempt
                Dictionary<int, Modification> tempModificationsDictionary = new(newModificationsDictionary);
                fillPosition = 0;
                // residueNums is a list containing array indices for each element of evaporatingBase
                // Once each amino acid is added, its index is removed from residueNums to prevent the same AA from being added 2x
                var residueNums = Enumerable.Range(0, evaporatingBase.Length).ToList();                
                characterCounter = 0;
                char[] tempNewBase = new char[newBase.Length];
                // Create a copy of the newBase character array for the scrambling attempt
                Array.Copy(newBase, tempNewBase, newBase.Length);

                // I am not sure why I need the second counter, but it always works when I have it
                int seqLength = this.BaseSequence.Length;
                while (fillPosition < seqLength && characterCounter < seqLength)
                {
                    residueNumsIndex = rand.Next(residueNums.Count);
                    extractPosition = residueNums[residueNumsIndex];
                    char targetAA = evaporatingBase[extractPosition];
                    residueNums.RemoveAt(residueNumsIndex);
                    if (targetAA != '0')
                    {
                        while (tempNewBase[fillPosition] != '0')
                        {
                            fillPosition++;
                        }
                        tempNewBase[fillPosition] = targetAA;
                        revisedAminoAcidOrder[fillPosition] = extractPosition;
                        if (this.AllModsOneIsNterminus.ContainsKey(extractPosition + 2))
                        {
                            tempModificationsDictionary.Add(fillPosition + 2, this.AllModsOneIsNterminus[extractPosition + 2]);
                        }
                        fillPosition++;
                    }
                    characterCounter ++;
                }
                scrambleAttempt++;
                /* 
                 * Any homology scoring mechanism can go here, percent identity is probably not the best
                 * In terms of generating a decoy sequence that will have a different mass spectrum than
                 * the original, it is far more important to vary the amino acids on the edges than 
                 * those in the middle. Changes on the edges will offset the entire b and y sequences
                 * leading to an effective decoy spectrum even if there is high identity in the middle of
                 * the sequence. Additionally, for peptides with a large amount of a certain amino acid,
                 * it will be very difficult to generate a low homology sequence.
                 */
                percentIdentity = GetPercentIdentity(tempNewBase, evaporatingBase, tempModificationsDictionary, this.AllModsOneIsNterminus);
                // Check that the percent identity is below the maximum identity threshold and set actual values to the temporary values
                if (percentIdentity < maxIdentity)
                {
                    newBase = tempNewBase;
                    newModificationsDictionary = tempModificationsDictionary;
                    // Code checking similarity between theoretical spectra could go here
                }

                // If max scrambles are reached, make the new sequence identical to the original to trigger mirroring
                else if (scrambleAttempt == maxScrambles)
                {
                    for(int j = 0; j < newBase.Length; j++)
                    {
                        if (newBase[j] == '0')
                        {
                            newBase[j] = evaporatingBase[j];
                        }
                    }
                }
            }
            

            string newBaseString = new string(newBase);

            var proteinSequence = this.Protein.BaseSequence;
            var aStringBuilder = new StringBuilder(proteinSequence);
            aStringBuilder.Remove(this.OneBasedStartResidueInProtein - 1, this.BaseSequence.Length);
            aStringBuilder.Insert(this.OneBasedStartResidueInProtein - 1, newBaseString);
            proteinSequence = aStringBuilder.ToString();

            Protein decoyProtein = new Protein(proteinSequence, "DECOY_" + this.Protein.Accession, null, new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), null, null, null, true);
            DigestionParams d = this.DigestionParams;
            // Creates a hash code corresponding to the target's sequence
            int targetHash = GetHashCode();
            PeptideWithSetModifications decoyPeptide;
            //Make the "peptideDescription" store the corresponding target's sequence
            if (newBaseString != this.BaseSequence)
            {
                decoyPeptide = new PeptideWithSetModifications(decoyProtein, d, this.OneBasedStartResidueInProtein, this.OneBasedEndResidueInProtein, this.CleavageSpecificityForFdrCategory, this.FullSequence, this.MissedCleavages, newModificationsDictionary, this.NumFixedMods, newBaseString);
                // Sets PairedTargetDecoyHash of the original target peptie to the hash hode of the decoy sequence           
                PairedTargetDecoyHash = decoyPeptide.GetHashCode();
                // Sets PairedTargetDecoyHash of the decoy peptide to the hash code of the target sequence
                decoyPeptide.PairedTargetDecoyHash = targetHash;
                return decoyPeptide;

            }
            else
            {
                //The reverse decoy procedure failed to create a PeptideWithSetModificatons with a different sequence. Therefore,
                //we retrun the mirror image peptide.
                decoyPeptide = this.GetPeptideMirror(revisedAminoAcidOrder);
                PairedTargetDecoyHash = decoyPeptide.GetHashCode();
                decoyPeptide.PairedTargetDecoyHash = targetHash;
                return decoyPeptide;
            }
        }
        
        /// <summary>
        /// Method to get the percent identity between two peptide sequences stored as char[]
        /// </summary>
        /// <param name="scrambledSequence">Character array of the scrambled sequence</param>
        /// <param name="unscrambledSequence">Character array of the unscrambled sequence</param> 
        /// <param name="scrambledMods">Dictionary containing the scrambled sequence's modifications</param>
        /// <param name="unscrambledMods">Dictionary containing the unscrambled sequence's modifications</param>
        /// <returns></returns>
        private static double GetPercentIdentity(char[] scrambledSequence, char[] unscrambledSequence, Dictionary<int, Modification> scrambledMods, Dictionary<int, Modification> unscrambledMods)
        {
            double rawScore = 0;
            int seqLength = scrambledSequence.Length;
            for(int i = 0; i < seqLength; i++)
            {
                if (scrambledSequence[i] == unscrambledSequence[i] || unscrambledSequence[i] == '0')
                {
                    Modification scrambledMod;
                    if (scrambledMods.TryGetValue(i + 2, out scrambledMod) && unscrambledSequence[i] != '0')
                    {
                        Modification unscrambledMod;
                        if (unscrambledMods.TryGetValue(i + 2, out unscrambledMod))
                        {
                            if (scrambledMod == unscrambledMod)
                            {
                                rawScore += 1;
                            }
                        }
                    }
                    else
                    {
                        rawScore += 1; 
                    }
                    
                }
            }
            return rawScore / seqLength;
        }
        
        //Returns a PeptideWithSetModifications mirror image. Used when reverse decoy sequence is same as target sequence
        public PeptideWithSetModifications GetPeptideMirror(int[] revisedOrderNisOne)
        {
            Dictionary<int, Modification> newModificationsDictionary = new Dictionary<int, Modification>();
            //Copy N-terminal modifications from target dictionary to decoy dictionary.
            if (this.AllModsOneIsNterminus.ContainsKey(1))
            {
                newModificationsDictionary.Add(1, this.AllModsOneIsNterminus[1]);
            }

            //First step is to reverse the position of all modifications except the mod on the peptide N-terminus.
            if (this.AllModsOneIsNterminus.Any())
            {
                foreach (var kvp in this.AllModsOneIsNterminus.Where(p => p.Key != 1).ToList())
                {
                    newModificationsDictionary.Add(this.BaseSequence.Length - kvp.Key + 3, kvp.Value);
                }
            }

            //Second step is to reverse the sequence.
            string newBaseString = new string(this.BaseSequence.Reverse().ToArray());

            var proteinSequence = this.Protein.BaseSequence;
            var aStringBuilder = new StringBuilder(proteinSequence);
            aStringBuilder.Remove(this.OneBasedStartResidueInProtein - 1, this.BaseSequence.Length);
            aStringBuilder.Insert(this.OneBasedStartResidueInProtein - 1, newBaseString);
            proteinSequence = aStringBuilder.ToString();

            Protein decoyProtein = new Protein(proteinSequence, "DECOY_" + this.Protein.Accession, null, new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), null, null, null, true);

            DigestionParams d = this.DigestionParams;

            //now fill in the revised amino acid order
            int oldStringPosition = this.BaseSequence.Length - 1;
            for (int i = 0; i < newBaseString.Length; i++)
            {
                revisedOrderNisOne[i] = oldStringPosition;
                oldStringPosition--;
            }

            //Make the "peptideDescription" store the corresponding target's sequence
            return new PeptideWithSetModifications(decoyProtein, d, this.OneBasedStartResidueInProtein, this.OneBasedEndResidueInProtein, this.CleavageSpecificityForFdrCategory, this.FullSequence, this.MissedCleavages, newModificationsDictionary, this.NumFixedMods, newBaseString);
        }
    }
}
