using System;
using MzLibUtil;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// Free-text instrument names, as a Thermo RAW file records them, resolved to PSI-MS instrument model terms.
    ///
    /// A name lookup alone is not enough. <see cref="ControlledVocabulary.TryGetByName"/> returns any term with that
    /// name or exact synonym, and the name a RAW file falls back to when it records no model is its device family,
    /// typically "MS" (ThermoRawFileReader.BuildInstrumentModel). psi-ms.obo imports UO, whose UO:0000028 "millisecond"
    /// has the exact synonym "ms", and the name index ignores case, so "MS" resolves to a unit of time. Written into an mzML's instrument configuration or an SDRF's
    /// comment[instrument], it would name an instrument nobody read. So a resolved term counts only when it descends
    /// from MS:1000031 "instrument model" and is not an abstract vendor branch ("Thermo Fisher Scientific instrument
    /// model"): the same branches the mzML reader refuses to report as a model (Mzml.FirstInstrumentModel).
    ///
    /// Shared by MzmlMethods and SdrfBuilder, so the two writers cannot disagree about what a RAW file's instrument is.
    ///
    /// A RAW file from the Q Exactive series names the model with the analyzer after it: "Q Exactive Orbitrap",
    /// "Q Exactive Plus Orbitrap", "Q Exactive HF Orbitrap", "Q Exactive HF-X Orbitrap". PSI-MS names the model alone
    /// ("Q Exactive HF", MS:1002523) and holds no synonym with the suffix, so the exact lookup fails, and every
    /// calibrated copy of such a file was written with the bare parent term. So a name that fails is tried once more
    /// without a trailing " Orbitrap", and the shorter name counts only when PSI-MS classes it as an Orbitrap
    /// instrument: the suffix claims an Orbitrap, and the ontology must agree ("LTQ Velos Orbitrap" does not become the
    /// ion trap "LTQ Velos").
    /// </summary>
    internal static class InstrumentModelCvTerms
    {
        /// <summary>MS:1000031 "instrument model", the parent of every instrument model term.</summary>
        internal const string InstrumentModelAccession = "MS:1000031";

        /// <summary>What a Q Exactive-series RAW file puts after the model's PSI-MS name.</summary>
        private const string OrbitrapSuffix = " Orbitrap";

        /// <summary>
        /// The PSI-MS instrument classes (children of MS:1003761 "instrument class") with an Orbitrap analyzer. A
        /// name that loses <see cref="OrbitrapSuffix"/> must resolve to a model in one of them.
        /// </summary>
        private static readonly string[] OrbitrapInstrumentClasses =
        {
            "MS:1003953", // orbitrap instrument (Exactive, Exactive Plus)
            "MS:1003768", // ion trap orbitrap instrument
            "MS:1003769", // quadrupole orbitrap instrument (the Q Exactive series, the Exploris series)
            "MS:1003770", // quadrupole ion trap orbitrap instrument
            "MS:1003771", // quadrupole orbitrap astral instrument
        };

        /// <summary>
        /// The PSI-MS instrument model with this name or exact synonym, ignoring case and whitespace, or with this name
        /// less a trailing " Orbitrap" when that names an Orbitrap instrument (see the class remarks). False, with a
        /// null <paramref name="model"/>, for a blank name, a name PSI-MS does not hold, and a name belonging to a
        /// term that is not a specific instrument model.
        /// </summary>
        internal static bool TryGetByName(string name, out CvParam model)
        {
            if (TryGetModel(name, out model))
                return true;

            string trimmed = name?.Trim();
            if (trimmed is null
                || !trimmed.EndsWith(OrbitrapSuffix, StringComparison.OrdinalIgnoreCase)
                || !TryGetModel(trimmed.Substring(0, trimmed.Length - OrbitrapSuffix.Length), out var shorter)
                || !Array.Exists(OrbitrapInstrumentClasses, c => ControlledVocabulary.PsiMs.IsA(shorter.Accession, c)))
                return false;

            model = shorter;
            return true;
        }

        /// <summary>The exact lookup: a name or exact synonym of a specific instrument model.</summary>
        private static bool TryGetModel(string name, out CvParam model)
        {
            model = null;
            if (!ControlledVocabulary.PsiMs.TryGetByName(name, out var term)
                || !ControlledVocabulary.PsiMs.IsA(term.Accession, InstrumentModelAccession)
                || term.Name.EndsWith("instrument model", StringComparison.OrdinalIgnoreCase))
                return false;

            model = term;
            return true;
        }
    }
}
