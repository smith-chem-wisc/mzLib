using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests
{
    /// <summary>
    /// Tests for the pinned, embedded PSI-MS and PRIDE controlled vocabularies.
    ///
    /// The point of pinning is that these answers do not move under the caller. Several assertions
    /// below therefore name specific accessions on purpose: if a snapshot update changes what
    /// "Q Exactive" or "HCD" resolves to, that is exactly the kind of change that should surface as
    /// a failing test and a deliberate decision, not as a silent shift in every file written
    /// afterwards.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestControlledVocabulary
    {
        [Test]
        public void PsiMs_LoadsAndReportsItsVersion()
        {
            var cv = ControlledVocabulary.PsiMs;

            Assert.That(cv.Label, Is.EqualTo("MS"));
            Assert.That(cv.Terms.Count, Is.GreaterThan(4000), "psi-ms.obo carries ~4,100 terms");
            Assert.That(cv.Version, Is.Not.Empty, "the data-version stamp is what makes the pin auditable");
        }

        [Test]
        public void Pride_LoadsAndReportsItsVersion()
        {
            var cv = ControlledVocabulary.Pride;

            Assert.That(cv.Label, Is.EqualTo("PRIDE"));
            Assert.That(cv.Terms.Count, Is.GreaterThan(900), "pride_cv.obo carries ~981 terms");
            Assert.That(cv.Version, Is.Not.Empty);
        }

        [Test]
        [TestCase("MS:1001911", "Q Exactive")]
        [TestCase("MS:1002732", "Orbitrap Fusion Lumos")]
        [TestCase("MS:1001251", "Trypsin")]
        [TestCase("MS:1000422", "beam-type collision-induced dissociation")]
        public void PsiMs_ResolvesByAccession(string accession, string expectedName)
        {
            Assert.That(ControlledVocabulary.PsiMs.TryGetByAccession(accession, out var term), Is.True);
            Assert.That(term.Name, Is.EqualTo(expectedName));
            Assert.That(term.CvLabel, Is.EqualTo("MS"));
            Assert.That(term.Value, Is.Empty, "a vocabulary term is a definition, not a measurement");
        }

        [Test]
        public void Pride_ResolvesTheAcquisitionMethodTerm()
        {
            // The corpus writes this from PRIDE, not PSI-MS -- which is why the smaller vocabulary
            // is not optional.
            Assert.That(ControlledVocabulary.Pride.TryGetByAccession("PRIDE:0000627", out var term), Is.True);
            Assert.That(term.Name, Is.EqualTo("Data-dependent acquisition"));
        }

        [Test]
        [TestCase("Q Exactive", "MS:1001911")]
        [TestCase("q exactive", "MS:1001911")]
        [TestCase("  Orbitrap   Fusion Lumos ", "MS:1002732")]
        public void PsiMs_ResolvesByNameCaseAndWhitespaceInsensitively(string name, string expected)
        {
            // This is the path that turns a Thermo RAW's free-text model into an accession.
            Assert.That(ControlledVocabulary.PsiMs.TryGetByName(name, out var term), Is.True);
            Assert.That(term.Accession, Is.EqualTo(expected));
        }

        [Test]
        public void UnknownAccessionAndNameResolveToFalse()
        {
            Assert.That(ControlledVocabulary.PsiMs.TryGetByAccession("MS:9999999", out _), Is.False);
            Assert.That(ControlledVocabulary.PsiMs.TryGetByName("not a real instrument", out _), Is.False);
            Assert.That(ControlledVocabulary.PsiMs.TryGetByAccession(null, out _), Is.False);
            Assert.That(ControlledVocabulary.PsiMs.TryGetByName("   ", out _), Is.False);
        }

        /// <summary>
        /// The obsolete filter is load-bearing, and this test now actually exercises it.
        ///
        /// The previous version asserted things about MS:1001911 "Q Exactive", which is not obsolete
        /// and has no obsolete namesake — so deleting the entire obsolete branch from the loader
        /// left it green. Sixteen obsolete names in psi-ms.obo shadow a live term and fifteen of
        /// them appear FIRST in file order, so without the filter "first writer wins" hands the
        /// lookup to the dead term every time.
        /// </summary>
        [Test]
        [TestCase("minute", "UO:0000031", "MS:1000038")]
        [TestCase("second", "UO:0000010", "MS:1000039")]
        [TestCase("dalton", "UO:0000221", "MS:1000212")]
        [TestCase("surface ionization", "MS:1000406", "MS:1000280")]
        public void ObsoleteTermsNeverWinANameLookup(string name, string liveAccession, string obsoleteAccession)
        {
            var cv = ControlledVocabulary.PsiMs;

            Assert.That(cv.TryGetByName(name, out var term), Is.True);
            Assert.That(term.Accession, Is.EqualTo(liveAccession),
                $"'{name}' resolved to the obsolete {obsoleteAccession}, which appears earlier in the file");

            // But the obsolete term must still be reachable BY ACCESSION: an old document
            // legitimately contains it and has to stay readable.
            Assert.That(cv.TryGetByAccession(obsoleteAccession, out var dead), Is.True);
            Assert.That(dead.Accession, Is.EqualTo(obsoleteAccession));
        }

        // ---------------- dissociation types ----------------

        [Test]
        [TestCase(DissociationType.HCD, "MS:1000422")]
        [TestCase(DissociationType.CID, "MS:1000133")]
        [TestCase(DissociationType.ETD, "MS:1000598")]
        [TestCase(DissociationType.ECD, "MS:1000250")]
        [TestCase(DissociationType.EThcD, "MS:1002631")]
        public void DissociationType_ResolvesToTheCanonicalTerm(DissociationType type, string expected)
        {
            Assert.That(DissociationTypeCvTerms.TryGetTerm(type, out var term), Is.True);
            Assert.That(term.Accession, Is.EqualTo(expected));
            Assert.That(term.Name, Is.Not.Empty, "the name comes from the pinned vocabulary, not a literal");
        }

        [Test]
        [TestCase("MS:1002481")]      // higher energy beam-type CID
        [TestCase("PRIDE:0000590")]   // PRIDE's HCD term
        [TestCase("MS:1000422")]      // the one we write
        public void Hcd_IsAcceptedUnderEverySpellingTheCorpusUses(string accession)
        {
            // Measured over the 950 corpus files carrying comment[dissociation method]: MS:1000422
            // in 273 files, MS:1002481 in 47, PRIDE:0000590 in 23. We write one and read all three,
            // so their files still join to ours.
            Assert.That(DissociationTypeCvTerms.TryGetDissociationType(accession, out var type), Is.True);
            Assert.That(type, Is.EqualTo(DissociationType.HCD));
        }

        /// <summary>
        /// Only four values genuinely have no PSI-MS term. An earlier revision claimed seven and
        /// was wrong about three of them: UVPD is MS:1003246 and NETD is MS:1003247 — both real,
        /// both used in the curated corpus, so refusing to write them was a joinability gap of my
        /// own making — and LowCID has MS:1002472 (trap-type CID), the structural counterpart to
        /// HCD's beam-type term. Asserting "no term exists" is a claim about the ontology, and it
        /// has to be checked against the ontology rather than assumed.
        /// </summary>
        [Test]
        [TestCase(DissociationType.aEPD)]
        [TestCase(DissociationType.Custom)]
        [TestCase(DissociationType.Autodetect)]
        [TestCase(DissociationType.AnyActivationType)]
        public void UnmappedDissociationTypes_ReturnFalseRatherThanAPlausibleNeighbour(DissociationType type)
        {
            // Custom/Autodetect/AnyActivationType are mzLib control values that name no method at
            // all. aEPD has no term: MS:1003294 "electron activated dissociation" is SCIEX
            // electron-beam EAD, a different technique. Returning a near-miss would be worse than
            // returning nothing, because a wrong accession silently joins to the wrong population.
            Assert.That(DissociationTypeCvTerms.TryGetTerm(type, out var term), Is.False);
            Assert.That(term, Is.Null);
        }

        /// <summary>
        /// Every mapped type, pinned to the NAME its accession carries — not merely to the fact
        /// that the accession resolves.
        ///
        /// Checking only that a term resolves is not enough, and this test exists because that gap
        /// let a real error through: SID was mapped to MS:1000406, which resolves perfectly well
        /// and is "surface ionization" (is_a MS:1000008 ionization type), not surface-induced
        /// dissociation (MS:1000136). A wrong-but-existing accession is the exact failure this
        /// class's own documentation warns about, and only asserting the name catches it.
        /// </summary>
        [Test]
        [TestCase(DissociationType.CID, "collision-induced dissociation")]
        [TestCase(DissociationType.HCD, "beam-type collision-induced dissociation")]
        [TestCase(DissociationType.ETD, "electron transfer dissociation")]
        [TestCase(DissociationType.ECD, "electron capture dissociation")]
        [TestCase(DissociationType.EThcD, "electron-transfer/higher-energy collision dissociation")]
        [TestCase(DissociationType.IRMPD, "infrared multiphoton dissociation")]
        [TestCase(DissociationType.PQD, "pulsed q dissociation")]
        [TestCase(DissociationType.BIRD, "blackbody infrared radiative dissociation")]
        [TestCase(DissociationType.SID, "surface-induced dissociation")]
        [TestCase(DissociationType.PD, "plasma desorption")]
        [TestCase(DissociationType.MPD, "photodissociation")]
        [TestCase(DissociationType.UVPD, "ultraviolet photodissociation")]
        [TestCase(DissociationType.NETD, "negative electron transfer dissociation")]
        [TestCase(DissociationType.LowCID, "trap-type collision-induced dissociation")]
        [TestCase(DissociationType.SORI, "sustained off-resonance irradiation")]
        [TestCase(DissociationType.PSD, "post-source decay")]
        [TestCase(DissociationType.ISCID, "in-source collision-induced dissociation")]
        public void EveryMappedTypeNamesTheRightConcept(DissociationType type, string expectedName)
        {
            Assert.That(DissociationTypeCvTerms.TryGetTerm(type, out var term), Is.True,
                $"{type} is listed as mapped but did not resolve");
            Assert.That(term.Name, Is.EqualTo(expectedName).IgnoreCase);
            Assert.That(term.Accession, Does.StartWith("MS:"));
        }

        [Test]
        public void TheMappedSetAndTheNamePinsAreInStep()
        {
            // If someone adds a type to the table without adding a name assertion above, the new
            // entry would go unverified — which is how the SID error survived in the first place.
            Assert.That(DissociationTypeCvTerms.Mapped.Count, Is.EqualTo(17),
                "add a TestCase to EveryMappedTypeNamesTheRightConcept for the new entry, then " +
                "update this count");
        }

        /// <summary>
        /// AGREEMENT WITH mzLib's OWN mzML READER. This is the test that would have caught the
        /// errors the name-pinning test could not.
        ///
        /// Pinning an accession to the name the OBO gives it is CIRCULAR: the ontology guarantees
        /// that link, so the assertion can only fail if the accession does not exist or the name
        /// was mistyped. It cannot detect the failure that actually happened three times here — a
        /// real, resolvable accession that names the wrong CONCEPT. SID was mapped to "surface
        /// ionization", PD to "photodissociation" (mzLib's PD is plasma desorption), and two read
        /// aliases pointed at "collision energy" and "Alkylation reagent". Every one of them
        /// resolved cleanly and passed a name-pinning test.
        ///
        /// Readers.Mzml.DissociationDictionary already encoded the right answers, in the same
        /// assembly, and disagreeing with it means an mzML this library reads cannot be described
        /// by the terms this library writes.
        /// </summary>
        [Test]
        public void AgreesWithTheMzmlReadersOwnAccessionTable()
        {
            var disagreements = new List<string>();

            foreach (var pair in Readers.Mzml.DissociationDictionary)
            {
                // MS:1000044 is the generic parent "dissociation method". Mzml maps it to Unknown
                // (Mzml.cs:112) because a file declaring only the parent has not said which method
                // it used. There is nothing for this table to agree or disagree with: it correctly
                // refuses to name a specific type, exactly as it does for the bare MS:1000031
                // instrument parent elsewhere.
                if (pair.Value == DissociationType.Unknown)
                    continue;

                if (!DissociationTypeCvTerms.TryGetDissociationType(pair.Key, out var ours))
                {
                    disagreements.Add($"{pair.Key}: mzML reader says {pair.Value}, we cannot read it");
                    continue;
                }
                if (ours != pair.Value)
                    disagreements.Add($"{pair.Key}: mzML reader says {pair.Value}, we say {ours}");
            }

            Assert.That(disagreements, Is.Empty,
                "the mzML reader and the CV table disagree: " + string.Join("; ", disagreements));
        }

        /// <summary>
        /// Read-aliases get the same name-pinning discipline as the primary table. They did not
        /// before, and both of the aliases that existed were wrong.
        /// </summary>
        [Test]
        [TestCase("MS:1002481", DissociationType.HCD, "higher energy beam-type collision-induced dissociation")]
        [TestCase("PRIDE:0000590", DissociationType.HCD, "HCD")]
        [TestCase("PRIDE:0000591", DissociationType.CID, "CID")]
        [TestCase("PRIDE:0000589", DissociationType.EThcD, "EThcD")]
        [TestCase("MS:1000433", DissociationType.LowCID, "low-energy collision-induced dissociation")]
        public void EveryReadAliasNamesTheRightConcept(string accession, DissociationType expected, string expectedName)
        {
            Assert.That(DissociationTypeCvTerms.TryGetDissociationType(accession, out var type), Is.True);
            Assert.That(type, Is.EqualTo(expected));

            // And confirm the accession really is the term we think, in whichever vocabulary owns it.
            var cv = accession.StartsWith("PRIDE:") ? ControlledVocabulary.Pride : ControlledVocabulary.PsiMs;
            Assert.That(cv.TryGetByAccession(accession, out var term), Is.True);
            Assert.That(term.Name, Is.EqualTo(expectedName).IgnoreCase);
        }

        [Test]
        public void CvLabelComesFromTheAccessionNotTheFile()
        {
            // psi-ms.obo imports 128 foreign terms and pride_cv.obo 19. Stamping the file's label on
            // them would emit cvRef="MS" for a UO term.
            Assert.That(ControlledVocabulary.PsiMs.TryGetByAccession("UO:0000266", out var electronvolt), Is.True);
            Assert.That(electronvolt.CvLabel, Is.EqualTo("UO"), "a unit term is not an MS term");

            Assert.That(ControlledVocabulary.Pride.TryGetByAccession("MS:1000044", out var method), Is.True);
            Assert.That(method.CvLabel, Is.EqualTo("MS"), "pride_cv.obo imports this MS term");
        }

        [Test]
        public void RoundTrip_TermToTypeAndBack()
        {
            foreach (var type in DissociationTypeCvTerms.Mapped)
            {
                Assert.That(DissociationTypeCvTerms.TryGetTerm(type, out var term), Is.True);
                Assert.That(DissociationTypeCvTerms.TryGetDissociationType(term.Accession, out var back), Is.True);

                Assert.That(back, Is.EqualTo(type),
                    $"{type} -> {term.Accession} -> {back}: every accession is now distinct, so a " +
                    "round trip must return the type it started from");
            }
        }
    }
}
