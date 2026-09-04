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
        /// Now case-insensitive, where the copy this replaced matched exactly. More permissive rather
        /// than different: every spelling the old switch accepted still resolves to the same value.
        /// </remarks>
        public static DissociationType? ModDissociationType(string modType) =>
            ModificationLoader.ParseDissociationType(modType);

        /// <inheritdoc cref="ModificationLoader.DiagnosticIonsAndNeutralLosses(string, Dictionary{DissociationType, List{double}})"/>
        public static Dictionary<DissociationType, List<double>> DiagnosticIonsAndNeutralLosses(string oneEntry,
            Dictionary<DissociationType, List<double>> dAndNDictionary) =>
            ModificationLoader.DiagnosticIonsAndNeutralLosses(oneEntry, dAndNDictionary);
    }
}
