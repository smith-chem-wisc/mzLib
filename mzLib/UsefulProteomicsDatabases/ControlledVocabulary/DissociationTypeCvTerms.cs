using System;
using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;
using MzLibUtil;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// Maps mzLib's <see cref="DissociationType"/> onto PSI-MS terms, and back.
    ///
    /// This table did not exist anywhere. The accessions were written in XML doc comments on the
    /// enum — where no code can reach them — and the only executable mapping was
    /// <c>Readers.Mzml.DissociationDictionary</c>, which is private to mzML reading, points the
    /// other way, and covers 15 of the enum's values.
    ///
    /// FOUR VALUES HAVE NO PSI-MS TERM and are absent here rather than being given a
    /// plausible-looking neighbour: aEPD, Custom, Autodetect and AnyActivationType. The last three
    /// are mzLib control values that do not name a method at all; aEPD genuinely has no term
    /// (MS:1003294 "electron activated dissociation" is SCIEX electron-beam EAD, a different
    /// thing). A caller that hits one gets false from <see cref="TryGetTerm"/> and should say so —
    /// writing a wrong accession is worse than writing none, because a wrong one silently joins to
    /// the wrong population.
    ///
    /// HCD is the one that needed deciding. Curated annotations spell it three ways, measured over
    /// the 950 corpus files carrying comment[dissociation method]: MS:1000422 in 273 files,
    /// MS:1002481 in 47, PRIDE:0000590 in 23. MS:1000422 ("beam-type collision-induced
    /// dissociation") is both the plurality and the term the specification's own examples use, so
    /// it is what we write; the others are accepted on read via
    /// <see cref="TryGetDissociationType"/> so documents using them still join to ours.
    ///
    /// EThcD is worth knowing about: MS:1002631 does NOT descend from MS:1000044, but from
    /// MS:1003181 "combined dissociation method". An ancestry-based validator built on MS:1000044
    /// will reject it, correctly-but-unhelpfully.
    /// </summary>
    public static class DissociationTypeCvTerms
    {
        private static readonly Dictionary<DissociationType, string> Accessions = new()
        {
            [DissociationType.CID] = "MS:1000133",     // collision-induced dissociation
            [DissociationType.HCD] = "MS:1000422",     // beam-type collision-induced dissociation
            [DissociationType.ETD] = "MS:1000598",     // electron transfer dissociation
            [DissociationType.ECD] = "MS:1000250",     // electron capture dissociation
            [DissociationType.EThcD] = "MS:1002631",   // Electron-Transfer/Higher-Energy Collision Dissociation
            [DissociationType.IRMPD] = "MS:1000262",   // infrared multiphoton dissociation
            [DissociationType.PQD] = "MS:1000599",     // pulsed q dissociation
            [DissociationType.BIRD] = "MS:1000242",    // blackbody infrared radiative dissociation
            [DissociationType.SID] = "MS:1000136",     // surface-induced dissociation
                                                       // NOT MS:1000406, which is "surface
                                                       // ionization" (is_a MS:1000008 ionization
                                                       // type) -- a real accession naming an
                                                       // entirely different concept.
            [DissociationType.PD] = "MS:1000134",      // plasma desorption -- NOT photodissociation.
                                                       // mzLib's PD means plasma desorption
                                                       // (DissociationType.cs:33, Mzml.cs:99); the
                                                       // photodissociation term belongs to MPD.
            [DissociationType.MPD] = "MS:1000435",     // photodissociation
            [DissociationType.SORI] = "MS:1000282",    // sustained off-resonance irradiation
            [DissociationType.PSD] = "MS:1000135",     // post-source decay
            [DissociationType.ISCID] = "MS:1001880",   // in-source collision-induced dissociation
            [DissociationType.UVPD] = "MS:1003246",    // ultraviolet photodissociation
            [DissociationType.NETD] = "MS:1003247",    // negative electron transfer dissociation
            [DissociationType.LowCID] = "MS:1002472",  // trap-type collision-induced dissociation --
                                                       // the structural counterpart to HCD's
                                                       // beam-type MS:1000422, and what mzLib's
                                                       // LowCID actually is (low-energy ion-trap
                                                       // CID; see DissociationTypeCollection).
        };

        /// <summary>
        /// Accessions accepted when READING that mean the same method as one we write. Reading
        /// generously and writing one canonical spelling is what lets a corpus stay joinable
        /// without inventing a private vocabulary.
        /// </summary>
        private static readonly Dictionary<string, DissociationType> ReadAliases = new(StringComparer.OrdinalIgnoreCase)
        {
            ["MS:1002481"] = DissociationType.HCD,       // higher energy beam-type CID
            ["PRIDE:0000590"] = DissociationType.HCD,    // PRIDE's own HCD term (is_a MS:1000044)
            ["PRIDE:0000591"] = DissociationType.CID,    // PRIDE's CID term, used by 2 corpus files
            ["PRIDE:0000589"] = DissociationType.EThcD,  // PRIDE's EThcD term, used by 2 corpus files
            ["MS:1000433"] = DissociationType.LowCID,    // low-energy CID, the other spelling
        };

        // PRIDE has NO ETD term -- its only dissociation children are EThcD (0000589), HCD (0000590),
        // CID (0000591) and ETciD (0000592). ETciD has no mzLib enum value, so a document using it
        // correctly gets false rather than being forced into a near-neighbour.

        // Two aliases were removed after an audit, and the reason is worth keeping:
        //   MS:1000045    "collision energy" -- is_a MS:1000510 precursor activation attribute. A
        //                 NUMERIC PARAMETER, not a method. It appears in essentially every mzML
        //                 precursor block, so aliasing it to CID would have misread HCD and ETD
        //                 runs as CID constantly.
        //   PRIDE:0000598 "Alkylation reagent" -- is_a PRIDE:0000026 Alkylation. Nothing to do with
        //                 ETD; it was simply wrong.
        // Both resolved cleanly against the vocabulary, which is exactly why "the accession exists"
        // is not evidence that a mapping is right.

        private static readonly Lazy<Dictionary<string, DissociationType>> ByAccession = new(() =>
        {
            var map = new Dictionary<string, DissociationType>(StringComparer.OrdinalIgnoreCase);
            foreach (var pair in Accessions)
                if (!map.ContainsKey(pair.Value))
                    map[pair.Value] = pair.Key;
            foreach (var alias in ReadAliases)
                map[alias.Key] = alias.Value;
            return map;
        });

        /// <summary>
        /// The PSI-MS term for a dissociation type, with the name taken from the pinned vocabulary
        /// rather than hardcoded, so it cannot drift from the accession. False when PSI-MS has no
        /// term for it — see the remarks on this class for which and why.
        /// </summary>
        public static bool TryGetTerm(DissociationType dissociationType, out CvParam term)
        {
            term = null;
            if (!Accessions.TryGetValue(dissociationType, out string accession))
                return false;

            if (ControlledVocabulary.PsiMs.TryGetByAccession(accession, out term))
                return true;

            // The accession is in this table but not in the embedded snapshot: the snapshot moved
            // and this table did not. Fail rather than invent a name.
            term = null;
            return false;
        }

        /// <summary>
        /// The dissociation type an accession names, accepting the read-aliases above.
        /// </summary>
        public static bool TryGetDissociationType(string accession, out DissociationType dissociationType)
        {
            dissociationType = DissociationType.Unknown;
            return !string.IsNullOrWhiteSpace(accession)
                   && ByAccession.Value.TryGetValue(accession.Trim(), out dissociationType);
        }

        /// <summary>
        /// The dissociation types this table can write a term for. Exposed so a caller can check
        /// coverage up front instead of discovering a gap mid-run.
        /// </summary>
        public static IReadOnlyCollection<DissociationType> Mapped => Accessions.Keys;
    }
}
