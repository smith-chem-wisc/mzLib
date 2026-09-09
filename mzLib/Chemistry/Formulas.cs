namespace Chemistry;

/// <summary>
/// Named chemical formulas shared across the solution. Single source of truth for the canonical
/// nucleic acid residue chemistry consumed by Transcriptomics.Nucleotide and the MODOMICS loader
/// (Omics.Modifications.IO.Modomics).
/// </summary>
public static class Formulas
{
    public static ChemicalFormula WaterChemicalFormula => ChemicalFormula.ParseFormula("H2O");
    public static ChemicalFormula CarbonDioxideChemicalFormula => ChemicalFormula.ParseFormula("CO2");

    /// <summary>
    /// Adenine in its N-glycosidic bonded form (free base less one H), i.e. the base as it is carried
    /// within an RNA or DNA residue.
    /// </summary>
    public static ChemicalFormula AdenineBaseChemicalFormula => ChemicalFormula.ParseFormula("C5H4N5");

    /// <summary>Cytosine, bonded form (free base less one H).</summary>
    public static ChemicalFormula CytosineBaseChemicalFormula => ChemicalFormula.ParseFormula("C4H4N3O1");

    /// <summary>Guanine, bonded form (free base less one H).</summary>
    public static ChemicalFormula GuanineBaseChemicalFormula => ChemicalFormula.ParseFormula("C5H4N5O1");

    /// <summary>Uracil, bonded form (free base less one H). Also the base carried by pseudouridine (Y).</summary>
    public static ChemicalFormula UracilBaseChemicalFormula => ChemicalFormula.ParseFormula("C4H3N2O2");

    /// <summary>Thymine, bonded form (free base less one H).</summary>
    public static ChemicalFormula ThymineBaseChemicalFormula => ChemicalFormula.ParseFormula("C5H5N2O2");

    /// <summary>
    /// Neutral adenosine (bonded adenine + ribose, condensed): the reference nucleoside that MODOMICS
    /// nucleoside-formula entries are stated against.
    /// </summary>
    public static ChemicalFormula AdenosineChemicalFormula => ChemicalFormula.ParseFormula("C10H13N5O4");

    /// <summary>Neutral cytidine (bonded cytosine + ribose, condensed).</summary>
    public static ChemicalFormula CytidineChemicalFormula => ChemicalFormula.ParseFormula("C9H13N3O5");

    /// <summary>Neutral guanosine (bonded guanine + ribose, condensed).</summary>
    public static ChemicalFormula GuanosineChemicalFormula => ChemicalFormula.ParseFormula("C10H13N5O5");

    /// <summary>Neutral uridine (bonded uracil + ribose, condensed).</summary>
    public static ChemicalFormula UridineChemicalFormula => ChemicalFormula.ParseFormula("C9H12N2O6");
}
