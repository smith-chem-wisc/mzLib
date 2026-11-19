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
    /// </summary>
    internal static class InstrumentModelCvTerms
    {
        /// <summary>MS:1000031 "instrument model", the parent of every instrument model term.</summary>
        internal const string InstrumentModelAccession = "MS:1000031";

        /// <summary>
        /// The PSI-MS instrument model with this name or exact synonym, ignoring case and whitespace. False, with a
        /// null <paramref name="model"/>, for a blank name, a name PSI-MS does not hold, and a name belonging to a
        /// term that is not a specific instrument model.
        /// </summary>
        internal static bool TryGetByName(string name, out CvParam model)
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
