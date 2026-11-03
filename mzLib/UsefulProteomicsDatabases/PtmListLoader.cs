using System;
using System.Collections.Generic;
using System.IO;
using MassSpectrometry;
using Omics.Modifications;
using Omics.Modifications.IO;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// Legacy wrapper for backwards compatibility.
    /// All functionality lives in <see cref="ModificationLoader"/>; every member here forwards to it.
    /// </summary>
    /// <remarks>
    /// The copy this replaced dropped BM and BL entries; they now load as
    /// <see cref="BackboneModification"/> and <see cref="BaseModification"/>.
    /// </remarks>
    [Obsolete("Use Omics.Modifications.ModificationLoader instead.")]
    public static class PtmListLoader
    {
        /// <inheritdoc cref="ModificationLoader.ReadModsFromFile(string, out List{ValueTuple{Modification, string}})"/>
        public static IEnumerable<Modification> ReadModsFromFile(string ptmListLocation,
            out List<(Modification, string)> filteredModificationsWithWarnings) =>
            ModificationLoader.ReadModsFromFile(ptmListLocation, out filteredModificationsWithWarnings);

        /// <inheritdoc cref="ModificationLoader.ReadModsFromFile(string, Dictionary{string, int}, out List{ValueTuple{Modification, string}})"/>
        public static IEnumerable<Modification> ReadModsFromFile(string ptmListLocation,
            Dictionary<string, int> formalChargesDictionary,
            out List<(Modification, string)> filteredModificationsWithWarnings) =>
            ModificationLoader.ReadModsFromFile(ptmListLocation, formalChargesDictionary,
                out filteredModificationsWithWarnings);

        /// <inheritdoc cref="ModificationLoader.ReadModsFromFile(StreamReader, Dictionary{string, int}, out List{ValueTuple{Modification, string}}, string)"/>
        public static IEnumerable<Modification> ReadModsFromFile(StreamReader uniprot_mods,
            Dictionary<string, int> formalChargesDictionary,
            out List<(Modification, string)> filteredModificationsWithWarnings, string? fileLocation = null) =>
            ModificationLoader.ReadModsFromFile(uniprot_mods, formalChargesDictionary,
                out filteredModificationsWithWarnings, fileLocation);

        /// <inheritdoc cref="ModificationLoader.ReadModsFromString(string, out List{ValueTuple{Modification, string}})"/>
        public static IEnumerable<Modification> ReadModsFromString(string storedModifications,
            out List<(Modification, string)> filteredModificationsWithWarnings) =>
            ModificationLoader.ReadModsFromString(storedModifications, out filteredModificationsWithWarnings);

        /// <inheritdoc cref="ModificationLoader.ParseDissociationType(string)"/>
        /// <remarks>
        /// Accepts more than the switch this replaced. That switch matched nine keys exactly (Any, CID,
        /// MPD, ECD, PQD, ETD, HCD, EThcD, Custom) and returned null for anything else; each of those
        /// still resolves to the same value. Matching is now case-insensitive, and every other
        /// <see cref="DissociationType"/> member name resolves to that member where it used to give null.
        /// </remarks>
        public static DissociationType? ModDissociationType(string modType) =>
            ModificationLoader.ParseDissociationType(modType);

        /// <inheritdoc cref="ModificationLoader.DiagnosticIonsAndNeutralLosses(string, Dictionary{DissociationType, List{double}})"/>
        public static Dictionary<DissociationType, List<double>> DiagnosticIonsAndNeutralLosses(string oneEntry,
            Dictionary<DissociationType, List<double>> dAndNDictionary) =>
            ModificationLoader.DiagnosticIonsAndNeutralLosses(oneEntry, dAndNDictionary);
    }
}
