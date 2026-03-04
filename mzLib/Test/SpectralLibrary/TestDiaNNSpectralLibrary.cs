// Test project setup:
//
//   dotnet add package NUnit --version 3.*
//   dotnet add package NUnit3TestAdapter
//   dotnet add package Microsoft.NET.Test.Sdk
//   dotnet add package Parquet.Net
//
// Target framework: net8.0
//
// These tests cover:
//   DiaNNModificationMapping     (Prompt 2)
//   DiaNNBinaryStructs           (Prompt 2)
//   DiaNNLibraryEntry            (Prompt 2)
//   DiaNNTsvSpectralLibrary      (Prompt 2)
//   DiaNNSpecLibReader           (Prompt 3)
//   DiaNNSpecLibWriter           (Prompt 4)
//   DiaNNSpecLibIndex            (Prompt 4)
//   DiaNNParquetSpectralLibrary  (Prompt 5)
//   Integration: TSV → Entry → Parquet → Entry round-trip
//   Integration: Entry → Binary → Entry round-trip
//   Integration: Index query correctness

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using NUnit.Framework;
using Omics.SpectrumMatch;
using Readers.SpectralLibrary.DiaNNSpectralLibrary;

namespace TestDiaNNSpectralLibrary
{
    // ═════════════════════════════════════════════════════════════════════════════════════════
    // Shared test data factory
    // ═════════════════════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// Produces reproducible synthetic library entries for use across all test classes.
    /// All entries are physically plausible (real-ish m/z values, normalized intensities).
    /// </summary>
    internal static class TestData
    {
        /// <summary>
        /// A minimal single-entry library: unmodified peptide PEPTIDER at charge 2.
        /// </summary>
        public static DiaNNLibraryEntry SimpleEntry() => new DiaNNLibraryEntry
        {
            ModifiedSequence = "_PEPTIDER_",
            StrippedSequence = "PEPTIDER",
            PrecursorCharge  = 2,
            PrecursorMz      = 481.2567,
            RetentionTime    = 42.3,
            IonMobility      = 0.0,
            ProteinId        = "P12345",
            ProteinName      = "TestProtein",
            GeneName         = "TPROT",
            IsProteotypic    = true,
            IsDecoy          = false,
            Fragments        = new List<DiaNNFragmentIon>
            {
                new DiaNNFragmentIon { Mz = 175.1190, Intensity = 1.0f,  IonType = 'y', SeriesNumber = 1, Charge = 1, LossType = "noloss" },
                new DiaNNFragmentIon { Mz = 303.1776, Intensity = 0.62f, IonType = 'y', SeriesNumber = 2, Charge = 1, LossType = "noloss" },
                new DiaNNFragmentIon { Mz = 400.2304, Intensity = 0.45f, IonType = 'y', SeriesNumber = 3, Charge = 1, LossType = "noloss" },
                new DiaNNFragmentIon { Mz = 129.1022, Intensity = 0.31f, IonType = 'b', SeriesNumber = 1, Charge = 1, LossType = "noloss" },
                new DiaNNFragmentIon { Mz = 226.1550, Intensity = 0.18f, IonType = 'b', SeriesNumber = 2, Charge = 1, LossType = "noloss" },
            }
        };

        /// <summary>
        /// An entry with oxidized methionine and a neutral loss fragment.
        /// </summary>
        public static DiaNNLibraryEntry OxidizedEntry() => new DiaNNLibraryEntry
        {
            ModifiedSequence = "_PEPTM[UniMod:35]IDER_",
            StrippedSequence = "PEPTMIDER",
            PrecursorCharge  = 2,
            PrecursorMz      = 537.2620,
            RetentionTime    = 38.7,
            IonMobility      = 0.0,
            ProteinId        = "P67890",
            ProteinName      = "OxidizedProtein",
            GeneName         = "OXPROT",
            IsProteotypic    = true,
            IsDecoy          = false,
            Fragments        = new List<DiaNNFragmentIon>
            {
                new DiaNNFragmentIon { Mz = 175.1190, Intensity = 1.0f,  IonType = 'y', SeriesNumber = 1, Charge = 1, LossType = "noloss" },
                new DiaNNFragmentIon { Mz = 303.1776, Intensity = 0.55f, IonType = 'y', SeriesNumber = 2, Charge = 1, LossType = "noloss" },
                new DiaNNFragmentIon { Mz = 450.2456, Intensity = 0.30f, IonType = 'y', SeriesNumber = 3, Charge = 1, LossType = "H2O"    },
            }
        };

        /// <summary>
        /// A decoy entry (reversed sequence, IsDecoy = true).
        /// </summary>
        public static DiaNNLibraryEntry DecoyEntry() => new DiaNNLibraryEntry
        {
            ModifiedSequence = "_REDITP[DECOY]_",
            StrippedSequence = "REDITP",
            PrecursorCharge  = 2,
            PrecursorMz      = 481.2567,
            RetentionTime    = 44.1,
            IonMobility      = 0.0,
            ProteinId        = "DECOY_P12345",
            ProteinName      = "DECOY_TestProtein",
            GeneName         = "DECOY_TPROT",
            IsProteotypic    = false,
            IsDecoy          = true,
            Fragments        = new List<DiaNNFragmentIon>
            {
                new DiaNNFragmentIon { Mz = 180.1025, Intensity = 1.0f,  IonType = 'y', SeriesNumber = 1, Charge = 1, LossType = "noloss" },
                new DiaNNFragmentIon { Mz = 277.1553, Intensity = 0.71f, IonType = 'y', SeriesNumber = 2, Charge = 1, LossType = "noloss" },
            }
        };

        /// <summary>
        /// A list of diverse entries suitable for library-level tests.
        /// Contains targets at different m/z values, charges, and RTs, plus one decoy.
        /// </summary>
        public static List<DiaNNLibraryEntry> SmallLibrary()
        {
            var entries = new List<DiaNNLibraryEntry>();

            // Spread across a representative m/z range (400–1200 Da)
            float[] mzValues = { 401.5f, 523.3f, 612.8f, 750.2f, 880.1f, 999.4f, 1150.6f };
            float[] rtValues = { 10.5f,  25.3f,  38.7f,  52.1f,  65.8f,  79.4f,  93.2f   };
            int[]   charges  = { 2,      2,      3,      2,      3,      2,      4       };

            for (int i = 0; i < mzValues.Length; i++)
            {
                entries.Add(new DiaNNLibraryEntry
                {
                    ModifiedSequence = $"_TESTPEPTIDE{i:D2}_",
                    StrippedSequence = $"TESTPEPTIDE{i:D2}",
                    PrecursorCharge  = charges[i],
                    PrecursorMz      = mzValues[i],
                    RetentionTime    = rtValues[i],
                    IonMobility      = 0.85f + i * 0.05f,
                    ProteinId        = $"PROT{i:D3}",
                    ProteinName      = $"Protein{i:D3}",
                    GeneName         = $"GENE{i:D3}",
                    IsProteotypic    = true,
                    IsDecoy          = false,
                    Fragments        = MakeSyntheticFragments(i)
                });
            }

            // Add one decoy at the end
            entries.Add(DecoyEntry());
            return entries;
        }

        private static List<DiaNNFragmentIon> MakeSyntheticFragments(int seed)
        {
            var frags = new List<DiaNNFragmentIon>();
            var rng   = new Random(seed);
            int count = rng.Next(4, 9);
            double baseMz = 150.0 + seed * 30.0;

            for (int i = 0; i < count; i++)
            {
                frags.Add(new DiaNNFragmentIon
                {
                    Mz           = (float)(baseMz + i * 101.05 + rng.NextDouble() * 5),
                    Intensity    = (float)(1.0 - i * 0.12 + rng.NextDouble() * 0.05),
                    IonType      = i % 2 == 0 ? 'y' : 'b',
                    SeriesNumber = i + 1,
                    Charge       = 1,
                    LossType     = "noloss"
                });
            }

            // Normalize so max = 1.0
            float maxInt = frags.Max(f => f.Intensity);
            foreach (var f in frags)
                f.Intensity /= maxInt;

            return frags;
        }

        /// <summary>A canonical DIA-NN TSV file for a two-precursor, minimal library.</summary>
        public static string MinimalTsvContent =>
            "ModifiedPeptide\tPrecursorCharge\tPrecursorMz\tTr_recalibrated\tProteinId\tFragmentMz\tRelativeIntensity\tFragmentType\tFragmentNumber\tFragmentCharge\tFragmentLossType\n" +
            "_PEPTIDER_\t2\t481.2567\t42.3\tP12345\t175.1190\t1.0000\ty\t1\t1\tnoloss\n" +
            "_PEPTIDER_\t2\t481.2567\t42.3\tP12345\t303.1776\t0.6200\ty\t2\t1\tnoloss\n" +
            "_PEPTM[UniMod:35]IDER_\t2\t537.2620\t38.7\tP67890\t175.1190\t1.0000\ty\t1\t1\tnoloss\n" +
            "_PEPTM[UniMod:35]IDER_\t2\t537.2620\t38.7\tP67890\t303.1776\t0.5500\ty\t2\t1\tnoloss\n";
    }

    // ═════════════════════════════════════════════════════════════════════════════════════════
    // 1. Modification mapping tests
    // ═════════════════════════════════════════════════════════════════════════════════════════

    [TestFixture]
    [Category("ModificationMapping")]
    public class DiaNNModificationMappingTests
    {
        [Test]
        public void DiaNNToMzLib_UnmodifiedSequence_RemovesUnderscoresOnly()
        {
            string result = DiaNNModificationMapping.DiaNNToMzLib("_PEPTIDER_");
            Assert.That(result, Is.EqualTo("PEPTIDER"));
        }

        [Test]
        public void DiaNNToMzLib_OxidizedMethionine_ConvertsCorrectly()
        {
            string result = DiaNNModificationMapping.DiaNNToMzLib("_PEPTM[UniMod:35]IDER_");
            Assert.That(result, Does.Contain("[Common Variable:Oxidation on M]"));
            Assert.That(result, Does.Not.Contain("_"));
            Assert.That(result, Does.Not.Contain("UniMod:35"));
        }

        [Test]
        public void DiaNNToMzLib_Carbamidomethyl_ConvertsCorrectly()
        {
            string result = DiaNNModificationMapping.DiaNNToMzLib("_PEPTC[UniMod:4]IDER_");
            Assert.That(result, Does.Contain("[Common Fixed:Carbamidomethyl on C]"));
            Assert.That(result, Does.Not.Contain("UniMod:4"));
        }

        [Test]
        public void MzLibToDiaNN_UnmodifiedSequence_AddsUnderscores()
        {
            string result = DiaNNModificationMapping.MzLibToDiaNN("PEPTIDER");
            Assert.That(result, Is.EqualTo("_PEPTIDER_"));
        }

        [Test]
        public void MzLibToDiaNN_OxidizedMethionine_ConvertsCorrectly()
        {
            string mzLib  = "PEPTM[Common Variable:Oxidation on M]IDER";
            string result = DiaNNModificationMapping.MzLibToDiaNN(mzLib);
            Assert.That(result, Does.Contain("[UniMod:35]"));
            Assert.That(result, Does.Not.Contain("Common Variable"));
        }

        [Test]
        public void RoundTrip_DiaNNToMzLibToDiaNN_PreservesSequence()
        {
            string original = "_PEPTM[UniMod:35]C[UniMod:4]IDER_";
            string mzLib    = DiaNNModificationMapping.DiaNNToMzLib(original);
            string restored = DiaNNModificationMapping.MzLibToDiaNN(mzLib);
            Assert.That(restored, Is.EqualTo(original));
        }

        [Test]
        public void GetStrippedSequence_RemovesAllModificationsAndUnderscores()
        {
            string result = DiaNNModificationMapping.GetStrippedSequence("_PEPTM[UniMod:35]C[UniMod:4]IDER_");
            Assert.That(result, Is.EqualTo("PEPTMCIDER"));
        }

        [Test]
        public void GetNeutralLossMass_KnownLosses_ReturnsCorrectMass()
        {
            Assert.That(DiaNNModificationMapping.GetNeutralLossMass("noloss"), Is.EqualTo(0.0).Within(1e-6));
            Assert.That(DiaNNModificationMapping.GetNeutralLossMass("H2O"),    Is.EqualTo(18.010565).Within(1e-5));
            Assert.That(DiaNNModificationMapping.GetNeutralLossMass("NH3"),    Is.EqualTo(17.026549).Within(1e-5));
        }

        [Test]
        public void GetNeutralLossName_KnownMasses_ReturnsCorrectName()
        {
            Assert.That(DiaNNModificationMapping.GetNeutralLossName(0.0),       Is.EqualTo("noloss"));
            Assert.That(DiaNNModificationMapping.GetNeutralLossName(18.010565), Is.EqualTo("H2O"));
        }
    }

    // ═════════════════════════════════════════════════════════════════════════════════════════
    // 2. Binary struct tests
    // ═════════════════════════════════════════════════════════════════════════════════════════

    [TestFixture]
    [Category("BinaryStructs")]
    public class DiaNNBinaryStructsTests
    {
        [Test]
        public void PrecursorRecord_SizeOf_IsPositiveAndAligned()
        {
            int size = DiaNNBinaryStructs.SizeOf.PrecursorRecord;
            Assert.That(size, Is.GreaterThan(0));
            // Minimum plausible: float(4) + short(2) + float(4)*2 + int(4)*3 + long(8) + short(2) + byte(2) + float(4) = 38
            Assert.That(size, Is.GreaterThanOrEqualTo(38));
            // Should be under 100 bytes — if it's larger something is wrong
            Assert.That(size, Is.LessThan(100));
        }

        [Test]
        public void FragmentRecord_SizeOf_IsExactly12Bytes()
        {
            // Per architecture doc: float(4) + float(4) + 4×byte(1) = 12
            int size = DiaNNBinaryStructs.SizeOf.FragmentRecord;
            Assert.That(size, Is.EqualTo(12));
        }

        [Test]
        public void IonTypeEncoding_RoundTrip_BAndY()
        {
            byte bByte = DiaNNBinaryStructs.IonTypeEncoding.ToByteFromChar('b');
            byte yByte = DiaNNBinaryStructs.IonTypeEncoding.ToByteFromChar('y');

            Assert.That(DiaNNBinaryStructs.IonTypeEncoding.ToCharFromByte(bByte), Is.EqualTo('b'));
            Assert.That(DiaNNBinaryStructs.IonTypeEncoding.ToCharFromByte(yByte), Is.EqualTo('y'));
        }

        [Test]
        public void LossTypeEncoding_RoundTrip_AllKnownLosses()
        {
            string[] known = { "noloss", "H2O", "NH3", "H3PO4", "HPO3" };
            foreach (var loss in known)
            {
                byte encoded = DiaNNBinaryStructs.LossTypeEncoding.ToByteFromString(loss);
                string decoded = DiaNNBinaryStructs.LossTypeEncoding.ToStringFromByte(encoded);
                Assert.That(decoded, Is.EqualTo(loss),
                    $"Round-trip failed for loss type '{loss}': got '{decoded}'");
            }
        }

        [Test]
        public void PrecursorRecord_IsSequentialPack1_NoPaddingBetweenFields()
        {
            // Verify Pack=1 is respected by checking that the struct has no hidden padding
            var layout = typeof(DiaNNBinaryStructs.PrecursorRecord)
                .GetCustomAttributes(typeof(StructLayoutAttribute), false)
                .Cast<StructLayoutAttribute>()
                .FirstOrDefault();

            Assert.That(layout, Is.Not.Null, "PrecursorRecord must have [StructLayout] attribute");
            Assert.That(layout!.Pack, Is.EqualTo(1), "PrecursorRecord must use Pack=1");
            Assert.That(layout.Value, Is.EqualTo(LayoutKind.Sequential));
        }
    }

    // ═════════════════════════════════════════════════════════════════════════════════════════
    // 3. LibraryEntry conversion tests
    // ═════════════════════════════════════════════════════════════════════════════════════════

    [TestFixture]
    [Category("LibraryEntry")]
    public class DiaNNLibraryEntryTests
    {
        [Test]
        public void ToLibrarySpectrum_SimpleEntry_HasCorrectPrecursorMz()
        {
            var entry    = TestData.SimpleEntry();
            var spectrum = entry.ToLibrarySpectrum();
            Assert.That(spectrum.PrecursorMz, Is.EqualTo(entry.PrecursorMz).Within(1e-4));
        }

        [Test]
        public void ToLibrarySpectrum_SimpleEntry_HasCorrectChargeState()
        {
            var entry    = TestData.SimpleEntry();
            var spectrum = entry.ToLibrarySpectrum();
            Assert.That(spectrum.ChargeState, Is.EqualTo(entry.PrecursorCharge));
        }

        [Test]
        public void ToLibrarySpectrum_SimpleEntry_HasCorrectFragmentCount()
        {
            var entry    = TestData.SimpleEntry();
            var spectrum = entry.ToLibrarySpectrum();
            Assert.That(spectrum.MatchedFragmentIons.Count, Is.EqualTo(entry.Fragments.Count));
        }

        [Test]
        public void ToLibrarySpectrum_AllFragmentMzValuesArePositive()
        {
            var entry    = TestData.SimpleEntry();
            var spectrum = entry.ToLibrarySpectrum();
            foreach (var ion in spectrum.MatchedFragmentIons)
                Assert.That(ion.Mz, Is.GreaterThan(0.0),
                    $"Fragment ion has non-positive m/z: {ion.Mz}");
        }

        [Test]
        public void ToLibrarySpectrum_DecoyEntry_IsDecoyFlagSet()
        {
            var entry    = TestData.DecoyEntry();
            var spectrum = entry.ToLibrarySpectrum();
            Assert.That(spectrum.IsDecoy, Is.True);
        }

        [Test]
        public void ToLibrarySpectrum_RetentionTimePreserved()
        {
            var entry    = TestData.SimpleEntry();
            var spectrum = entry.ToLibrarySpectrum();
            Assert.That(spectrum.RetentionTime, Is.EqualTo(entry.RetentionTime).Within(1e-3));
        }

        [Test]
        public void FromLibrarySpectrum_ThenToLibrarySpectrum_PrecursorMzPreserved()
        {
            var original  = TestData.SimpleEntry().ToLibrarySpectrum();
            var wrapped   = DiaNNLibraryEntry.FromLibrarySpectrum(original);
            var converted = wrapped.ToLibrarySpectrum();
            Assert.That(converted.PrecursorMz, Is.EqualTo(original.PrecursorMz).Within(1e-4));
        }

        [Test]
        public void FromLibrarySpectrum_ThenToLibrarySpectrum_FragmentCountPreserved()
        {
            var original  = TestData.SimpleEntry().ToLibrarySpectrum();
            var wrapped   = DiaNNLibraryEntry.FromLibrarySpectrum(original);
            var converted = wrapped.ToLibrarySpectrum();
            Assert.That(converted.MatchedFragmentIons.Count,
                Is.EqualTo(original.MatchedFragmentIons.Count));
        }

        [Test]
        public void Name_MatchesSequenceSlashCharge()
        {
            var entry    = TestData.SimpleEntry();
            var spectrum = entry.ToLibrarySpectrum();
            // mzLib LibrarySpectrum.Name = Sequence + "/" + ChargeState
            StringAssert.Contains("/2", spectrum.Name);
        }
    }

    // ═════════════════════════════════════════════════════════════════════════════════════════
    // 4. TSV reader/writer tests
    // ═════════════════════════════════════════════════════════════════════════════════════════

    [TestFixture]
    [Category("Tsv")]
    public class DiaNNTsvTests
    {
        private string _tempDir = null!;

        [SetUp]
        public void SetUp() => _tempDir = Path.Combine(Path.GetTempPath(), "DiaNNTsvTests_" + Guid.NewGuid());
        [TearDown]
        public void TearDown() { if (Directory.Exists(_tempDir)) Directory.Delete(_tempDir, recursive: true); }

        private string TempFile(string name) { Directory.CreateDirectory(_tempDir); return Path.Combine(_tempDir, name); }

        [Test]
        public void WriteThenRead_SimpleEntry_RoundTripsCorrectly()
        {
            var path    = TempFile("simple.tsv");
            var entries = new List<DiaNNLibraryEntry> { TestData.SimpleEntry() };

            DiaNNTsvSpectralLibrary.WriteTsv(path, entries);
            var loaded = DiaNNTsvSpectralLibrary.ReadTsv(path);

            Assert.That(loaded.Count, Is.EqualTo(1));
            Assert.That(loaded[0].PrecursorCharge, Is.EqualTo(2));
            Assert.That(loaded[0].PrecursorMz, Is.EqualTo(481.2567).Within(1e-3));
            Assert.That(loaded[0].Fragments.Count, Is.EqualTo(TestData.SimpleEntry().Fragments.Count));
        }

        [Test]
        public void WriteThenRead_MultipleEntries_AllRoundTrip()
        {
            var path    = TempFile("multi.tsv");
            var entries = TestData.SmallLibrary();

            DiaNNTsvSpectralLibrary.WriteTsv(path, entries);
            var loaded = DiaNNTsvSpectralLibrary.ReadTsv(path);

            Assert.That(loaded.Count, Is.EqualTo(entries.Count));
        }

        [Test]
        public void WriteThenRead_RetentionTimesPreserved()
        {
            var path    = TempFile("rt.tsv");
            var entries = TestData.SmallLibrary().Take(3).ToList();

            DiaNNTsvSpectralLibrary.WriteTsv(path, entries);
            var loaded = DiaNNTsvSpectralLibrary.ReadTsv(path);

            for (int i = 0; i < entries.Count; i++)
                Assert.That(loaded[i].RetentionTime,
                    Is.EqualTo(entries[i].RetentionTime).Within(1e-3),
                    $"RT mismatch at index {i}");
        }

        [Test]
        public void WriteThenRead_DecoyFlag_Preserved()
        {
            var path    = TempFile("decoy.tsv");
            var entries = new List<DiaNNLibraryEntry> { TestData.SimpleEntry(), TestData.DecoyEntry() };

            DiaNNTsvSpectralLibrary.WriteTsv(path, entries);
            var loaded = DiaNNTsvSpectralLibrary.ReadTsv(path);

            var target = loaded.Single(e => !e.IsDecoy);
            var decoy  = loaded.Single(e => e.IsDecoy);
            Assert.That(target.IsDecoy, Is.False);
            Assert.That(decoy.IsDecoy,  Is.True);
        }

        [Test]
        public void ReadMinimalTsvContent_ParsesTwoPrecursors()
        {
            var path = TempFile("minimal.tsv");
            File.WriteAllText(path, TestData.MinimalTsvContent);

            var loaded = DiaNNTsvSpectralLibrary.ReadTsv(path);

            Assert.That(loaded.Count, Is.EqualTo(2),
                "Expected exactly 2 precursors (PEPTIDER and the oxidized variant)");
        }

        [Test]
        public void ReadMinimalTsvContent_FirstPrecursorHasTwoFragments()
        {
            var path = TempFile("minimal2.tsv");
            File.WriteAllText(path, TestData.MinimalTsvContent);

            var loaded  = DiaNNTsvSpectralLibrary.ReadTsv(path);
            var peptider = loaded.First(e => e.ModifiedSequence.Contains("PEPTIDER") && !e.ModifiedSequence.Contains("UniMod"));

            Assert.That(peptider.Fragments.Count, Is.EqualTo(2));
        }

        [Test]
        public void ReadTsv_WithSpectronuautColumnAlias_ParsesSuccessfully()
        {
            // Use Spectronaut-style column name "ModifiedPeptideSequence" instead of "ModifiedPeptide"
            string tsvContent =
                "ModifiedPeptideSequence\tPrecursorCharge\tPrecursorMz\tNormalizedRetentionTime\tProteinId\tFragmentMz\tLibraryIntensity\tFragmentType\tFragmentNumber\tFragmentCharge\n" +
                "_PEPTIDER_\t2\t481.2567\t42.3\tP12345\t175.1190\t1.0\ty\t1\t1\n";

            var path = TempFile("spectronaut.tsv");
            File.WriteAllText(path, tsvContent);

            Assert.DoesNotThrow(() =>
            {
                var loaded = DiaNNTsvSpectralLibrary.ReadTsv(path);
                Assert.That(loaded.Count, Is.EqualTo(1));
            });
        }

        [Test]
        public void WriteFromLibrarySpectra_ThenRead_PrecursorCountMatches()
        {
            var spectra = TestData.SmallLibrary()
                .Take(3)
                .Select(e => e.ToLibrarySpectrum())
                .ToList();
            var path = TempFile("from_spectra.tsv");

            DiaNNTsvSpectralLibrary.WriteTsvFromLibrarySpectra(path, spectra);
            var loaded = DiaNNTsvSpectralLibrary.ReadTsv(path);

            Assert.That(loaded.Count, Is.EqualTo(spectra.Count));
        }

        [Test]
        public void ReadTsv_MissingRequiredColumn_ThrowsInvalidDataException()
        {
            // TSV without required FragmentMz column
            string badTsv =
                "ModifiedPeptide\tPrecursorCharge\tPrecursorMz\tTr_recalibrated\n" +
                "_PEPTIDER_\t2\t481.26\t42.3\n";

            var path = TempFile("bad.tsv");
            File.WriteAllText(path, badTsv);

            Assert.Throws<InvalidDataException>(() => DiaNNTsvSpectralLibrary.ReadTsv(path));
        }
    }

    // ═════════════════════════════════════════════════════════════════════════════════════════
    // 5. Binary writer + reader round-trip tests
    // ═════════════════════════════════════════════════════════════════════════════════════════

    [TestFixture]
    [Category("Binary")]
    public class DiaNNBinaryRoundTripTests
    {
        private string _tempDir = null!;

        [SetUp]
        public void SetUp() => _tempDir = Path.Combine(Path.GetTempPath(), "DiaNNBinaryTests_" + Guid.NewGuid());
        [TearDown]
        public void TearDown() { if (Directory.Exists(_tempDir)) Directory.Delete(_tempDir, recursive: true); }

        private string TempFile(string name) { Directory.CreateDirectory(_tempDir); return Path.Combine(_tempDir, name); }

        [Test]
        public void WriteRead_SingleEntry_PrecursorMzPreserved()
        {
            var path   = TempFile("single.speclib");
            var entry  = TestData.SimpleEntry();
            DiaNNSpecLibWriter.Write(path, new[] { entry });

            using var reader = new DiaNNSpecLibReader(path);
            var loaded = reader.ReadAllEntries().ToList();

            Assert.That(loaded.Count, Is.EqualTo(1));
            Assert.That(loaded[0].PrecursorMz, Is.EqualTo(entry.PrecursorMz).Within(1e-3));
        }

        [Test]
        public void WriteRead_SingleEntry_ChargePreserved()
        {
            var path  = TempFile("charge.speclib");
            var entry = TestData.SimpleEntry();
            DiaNNSpecLibWriter.Write(path, new[] { entry });

            using var reader = new DiaNNSpecLibReader(path);
            var loaded = reader.ReadAllEntries().ToList();

            Assert.That(loaded[0].PrecursorCharge, Is.EqualTo(entry.PrecursorCharge));
        }

        [Test]
        public void WriteRead_SingleEntry_FragmentCountPreserved()
        {
            var path  = TempFile("fragcount.speclib");
            var entry = TestData.SimpleEntry();
            DiaNNSpecLibWriter.Write(path, new[] { entry });

            using var reader = new DiaNNSpecLibReader(path);
            var loaded = reader.ReadAllEntries().ToList();

            Assert.That(loaded[0].Fragments.Count, Is.EqualTo(entry.Fragments.Count));
        }

        [Test]
        public void WriteRead_SingleEntry_FragmentMzValuesPreserved()
        {
            var path  = TempFile("fragmz.speclib");
            var entry = TestData.SimpleEntry();
            DiaNNSpecLibWriter.Write(path, new[] { entry });

            using var reader = new DiaNNSpecLibReader(path);
            var loaded = reader.ReadAllEntries().ToList();

            var origMzs   = entry.Fragments.Select(f => f.Mz).OrderBy(x => x).ToList();
            var loadedMzs = loaded[0].Fragments.Select(f => f.Mz).OrderBy(x => x).ToList();

            for (int i = 0; i < origMzs.Count; i++)
                Assert.That(loadedMzs[i], Is.EqualTo(origMzs[i]).Within(1e-3),
                    $"Fragment {i} m/z mismatch: expected {origMzs[i]:F4}, got {loadedMzs[i]:F4}");
        }

        [Test]
        public void WriteRead_SmallLibrary_AllEntriesRoundTrip()
        {
            var path    = TempFile("library.speclib");
            var entries = TestData.SmallLibrary();
            DiaNNSpecLibWriter.Write(path, entries);

            using var reader = new DiaNNSpecLibReader(path);
            var loaded = reader.ReadAllEntries().ToList();

            Assert.That(loaded.Count, Is.EqualTo(entries.Count));
        }

        [Test]
        public void WriteRead_DecoyFlag_Preserved()
        {
            var path    = TempFile("decoy.speclib");
            var entries = new List<DiaNNLibraryEntry> { TestData.SimpleEntry(), TestData.DecoyEntry() };
            DiaNNSpecLibWriter.Write(path, entries);

            using var reader = new DiaNNSpecLibReader(path);
            var loaded = reader.ReadAllEntries().ToList();

            Assert.That(loaded.Any(e => e.IsDecoy),  Is.True,  "Decoy entry missing after round-trip");
            Assert.That(loaded.Any(e => !e.IsDecoy), Is.True,  "Target entry missing after round-trip");
        }

        [Test]
        public void WriteRead_RetentionTimesPreserved_SmallLibrary()
        {
            var path    = TempFile("rt.speclib");
            var entries = TestData.SmallLibrary().Take(4).ToList();
            DiaNNSpecLibWriter.Write(path, entries);

            using var reader = new DiaNNSpecLibReader(path);
            var loaded = reader.ReadAllEntries()
                .OrderBy(e => e.PrecursorMz)
                .ToList();
            var original = entries.OrderBy(e => e.PrecursorMz).ToList();

            for (int i = 0; i < original.Count; i++)
                Assert.That(loaded[i].RetentionTime,
                    Is.EqualTo(original[i].RetentionTime).Within(1e-2),
                    $"RT mismatch at index {i}");
        }

        [Test]
        public void Write_EmptyList_ThrowsArgumentException()
        {
            var path = TempFile("empty.speclib");
            Assert.Throws<ArgumentException>(() =>
                DiaNNSpecLibWriter.Write(path, new List<DiaNNLibraryEntry>()));
        }

        [Test]
        public void ValidateEntries_ValidLibrary_ReturnsNoErrors()
        {
            var errors = DiaNNSpecLibWriter.ValidateEntries(TestData.SmallLibrary());
            Assert.That(errors, Is.Empty, string.Join("\n", errors));
        }

        [Test]
        public void ValidateEntries_DuplicateKey_ReportsError()
        {
            var duplicate = new List<DiaNNLibraryEntry>
            {
                TestData.SimpleEntry(),
                TestData.SimpleEntry() // exact duplicate
            };
            var errors = DiaNNSpecLibWriter.ValidateEntries(duplicate);
            Assert.That(errors.Any(e => e.Contains("Duplicate")), Is.True);
        }

        [Test]
        public void EstimateFileSize_IsReasonablyAccurate()
        {
            var path    = TempFile("estimate.speclib");
            var entries = TestData.SmallLibrary();

            long estimated = DiaNNSpecLibWriter.EstimateFileSize(entries);
            DiaNNSpecLibWriter.Write(path, entries);
            long actual = new FileInfo(path).Length;

            // Estimate should be within 2x of actual (strings may be overestimated without dedup)
            Assert.That(estimated, Is.LessThanOrEqualTo(actual * 2),
                $"Estimate {estimated:N0} is more than 2x actual {actual:N0}");
            Assert.That(estimated, Is.GreaterThan(0));
        }

        [Test]
        public void IsSpecLibFile_ValidFile_ReturnsTrue()
        {
            var path    = TempFile("valid.speclib");
            var entries = new List<DiaNNLibraryEntry> { TestData.SimpleEntry() };
            DiaNNSpecLibWriter.Write(path, entries);

            Assert.That(DiaNNSpecLibReader.IsSpecLibFile(path), Is.True);
        }

        [Test]
        public void IsSpecLibFile_TsvFile_ReturnsFalse()
        {
            var path = TempFile("not_speclib.tsv");
            File.WriteAllText(path, TestData.MinimalTsvContent);

            Assert.That(DiaNNSpecLibReader.IsSpecLibFile(path), Is.False);
        }

        [Test]
        public void WriteRead_FragmentsAreSortedByMz_AfterRoundTrip()
        {
            var path  = TempFile("sorted.speclib");
            var entry = TestData.SimpleEntry();
            // Deliberately write fragments in reverse m/z order
            entry.Fragments = entry.Fragments.OrderByDescending(f => f.Mz).ToList();
            DiaNNSpecLibWriter.Write(path, new[] { entry });

            using var reader = new DiaNNSpecLibReader(path);
            var loaded = reader.ReadAllEntries().ToList();
            var mzs    = loaded[0].Fragments.Select(f => f.Mz).ToList();

            for (int i = 1; i < mzs.Count; i++)
                Assert.That(mzs[i], Is.GreaterThanOrEqualTo(mzs[i - 1]),
                    $"Fragments not sorted at position {i}: {mzs[i - 1]:F4} > {mzs[i]:F4}");
        }
    }

    // ═════════════════════════════════════════════════════════════════════════════════════════
    // 6. Index tests
    // ═════════════════════════════════════════════════════════════════════════════════════════

    [TestFixture]
    [Category("Index")]
    public class DiaNNSpecLibIndexTests
    {
        private DiaNNSpecLibIndex _index = null!;
        private List<DiaNNLibraryEntry> _entries = null!;

        [SetUp]
        public void SetUp()
        {
            _entries = TestData.SmallLibrary();
            _index   = DiaNNSpecLibIndex.BuildFromEntries(_entries);
        }

        [TearDown]
        public void TearDown() => _index?.Dispose();

        // ── Range query tests ────────────────────────────────────────────────────────────

        [Test]
        public void GetPrecursorsInMzRange_CoverAll_ReturnsAllEntries()
        {
            var results = _index.GetPrecursorsInMzRange(0, 2000);
            Assert.That(results.Length, Is.EqualTo(_entries.Count));
        }

        [Test]
        public void GetPrecursorsInMzRange_EmptyRange_ReturnsEmpty()
        {
            var results = _index.GetPrecursorsInMzRange(0, 1);
            Assert.That(results.Length, Is.EqualTo(0));
        }

        [Test]
        public void GetPrecursorsInMzRange_SinglePrecursor_ReturnsCorrectEntry()
        {
            // First entry in SmallLibrary has m/z ≈ 401.5
            var results = _index.GetPrecursorsInMzRange(400.0, 410.0);
            Assert.That(results.Length, Is.EqualTo(1));
            Assert.That(results[0].PrecursorMz, Is.EqualTo(401.5f).Within(0.5f));
        }

        [Test]
        public void GetPrecursorsInMzRange_ResultsSortedByMz()
        {
            var results = _index.GetPrecursorsInMzRange(0, 2000);
            for (int i = 1; i < results.Length; i++)
                Assert.That(results[i].PrecursorMz, Is.GreaterThanOrEqualTo(results[i - 1].PrecursorMz),
                    $"Results not sorted at index {i}");
        }

        // ── Multi-dimensional window query tests ─────────────────────────────────────────

        [Test]
        public void GetPrecursorsInWindow_ChargeFilter_ExcludesWrongCharges()
        {
            using var scope = _index.GetPrecursorsInWindow(0, 2000, charge: 2);
            Assert.That(scope.Entries.ToArray().All(e => e.Charge == 2), Is.True);
        }

        [Test]
        public void GetPrecursorsInWindow_DecoyFilter_ExcludesDecoys()
        {
            using var scope = _index.GetPrecursorsInWindow(0, 2000, includeDecoys: false);
            Assert.That(scope.Entries.ToArray().Any(e => e.IsDecoy), Is.False);
        }

        [Test]
        public void GetPrecursorsInWindow_RtFilter_ConstrainsResults()
        {
            // RT range 20–60 should hit entries with RT 25.3, 38.7, and 52.1 from SmallLibrary
            using var scope = _index.GetPrecursorsInWindow(0, 2000, minRt: 20.0, maxRt: 60.0);
            Assert.That(scope.Count, Is.GreaterThanOrEqualTo(1));
            Assert.That(scope.Entries.ToArray().All(e => e.RetentionTime >= 20.0f && e.RetentionTime <= 60.0f),
                Is.True);
        }

        [Test]
        public void GetPrecursorsInWindow_EmptyResult_IsEmptyScope()
        {
            using var scope = _index.GetPrecursorsInWindow(9999, 10000);
            Assert.That(scope.IsEmpty, Is.True);
        }

        // ── DDA-style lookup ─────────────────────────────────────────────────────────────

        [Test]
        public void TryGetSpectrum_KnownEntry_ReturnsTrue()
        {
            var first = _entries[0];
            bool found = _index.TryGetSpectrum(first.ModifiedSequence!, first.PrecursorCharge, out var entry);

            Assert.That(found, Is.True);
            Assert.That(entry, Is.Not.Null);
        }

        [Test]
        public void TryGetSpectrum_UnknownSequence_ReturnsFalse()
        {
            bool found = _index.TryGetSpectrum("_DOESNOTEXIST_", 2, out _);
            Assert.That(found, Is.False);
        }

        [Test]
        public void TryGetLibrarySpectrum_KnownEntry_ReturnsLibrarySpectrum()
        {
            var first = _entries[0];
            bool found = _index.TryGetLibrarySpectrum(
                DiaNNModificationMapping.DiaNNToMzLib(first.ModifiedSequence!),
                first.PrecursorCharge,
                out var spectrum);

            Assert.That(found, Is.True);
            Assert.That(spectrum, Is.Not.Null);
            Assert.That(spectrum!.PrecursorMz, Is.EqualTo(first.PrecursorMz).Within(1e-3));
        }

        // ── Statistics ───────────────────────────────────────────────────────────────────

        [Test]
        public void GetStatistics_TotalPrecursorCount_MatchesInput()
        {
            var stats = _index.GetStatistics();
            Assert.That(stats.TotalPrecursors, Is.EqualTo(_entries.Count));
        }

        [Test]
        public void GetStatistics_DecoyCount_IsOne()
        {
            var stats = _index.GetStatistics();
            Assert.That(stats.DecoyPrecursors, Is.EqualTo(1));
        }

        [Test]
        public void GetStatistics_MzBoundsAreReasonable()
        {
            var stats = _index.GetStatistics();
            Assert.That(stats.MinPrecursorMz, Is.LessThan(stats.MaxPrecursorMz));
            Assert.That(stats.MinPrecursorMz, Is.GreaterThan(0));
        }

        // ── RT calibration ───────────────────────────────────────────────────────────────

        [Test]
        public void WithCalibratedRetentionTimes_IdentityTransform_LeavesRtUnchanged()
        {
            using var calibrated = _index.WithCalibratedRetentionTimes(slope: 1.0, intercept: 0.0);
            var original   = _index.AllIndexEntries.ToArray();
            var calEntries = calibrated.AllIndexEntries.ToArray();

            Assert.That(calEntries.Length, Is.EqualTo(original.Length));
            for (int i = 0; i < original.Length; i++)
                Assert.That(calEntries.Any(e => Math.Abs(e.RetentionTime - original[i].RetentionTime) < 0.01f),
                    Is.True,
                    $"RT {original[i].RetentionTime} not found after identity calibration");
        }

        [Test]
        public void WithCalibratedRetentionTimes_ScaledTransform_RtValuesScaled()
        {
            float firstRt = _index.AllIndexEntries[0].RetentionTime;
            using var calibrated = _index.WithCalibratedRetentionTimes(slope: 2.0, intercept: 5.0);

            // At least one RT should be approximately 2× + 5 of the original
            bool anyMatch = calibrated.AllIndexEntries.ToArray()
                .Any(e => Math.Abs(e.RetentionTime - (2.0f * firstRt + 5.0f)) < 0.1f);
            Assert.That(anyMatch, Is.True);
        }

        // ── QueryScope disposal ──────────────────────────────────────────────────────────

        [Test]
        public void QueryScope_Empty_IsEmptyAndSafeToDispose()
        {
            var scope = QueryScope.Empty;
            Assert.That(scope.IsEmpty, Is.True);
            Assert.DoesNotThrow(() => scope.Dispose());
        }

        [Test]
        public void QueryScope_AfterDispose_EntriesAreEmpty()
        {
            using var scope = _index.GetPrecursorsInWindow(0, 2000);
            int countBeforeDispose = scope.Count;
            Assert.That(countBeforeDispose, Is.GreaterThan(0));

            scope.Dispose();
            Assert.That(scope.Entries.Length, Is.EqualTo(0));
        }

        // ── Dispose safety ───────────────────────────────────────────────────────────────

        [Test]
        public void Dispose_CallingTwice_DoesNotThrow()
        {
            var index = DiaNNSpecLibIndex.BuildFromEntries(
                new[] { TestData.SimpleEntry() });
            index.Dispose();
            Assert.DoesNotThrow(() => index.Dispose());
        }

        [Test]
        public void AfterDispose_AnyQuery_ThrowsObjectDisposedException()
        {
            var index = DiaNNSpecLibIndex.BuildFromEntries(
                new[] { TestData.SimpleEntry() });
            index.Dispose();
            Assert.Throws<ObjectDisposedException>(() =>
                _ = index.GetPrecursorsInMzRange(0, 2000));
        }
    }

    // ═════════════════════════════════════════════════════════════════════════════════════════
    // 7. Parquet reader/writer tests
    // ═════════════════════════════════════════════════════════════════════════════════════════

    [TestFixture]
    [Category("Parquet")]
    public class DiaNNParquetTests
    {
        private string _tempDir = null!;

        [SetUp]
        public void SetUp() => _tempDir = Path.Combine(Path.GetTempPath(), "DiaNNParquetTests_" + Guid.NewGuid());
        [TearDown]
        public void TearDown() { if (Directory.Exists(_tempDir)) Directory.Delete(_tempDir, recursive: true); }

        private string TempFile(string name) { Directory.CreateDirectory(_tempDir); return Path.Combine(_tempDir, name); }

        [Test]
        public void WriteRead_SimpleEntry_PrecursorMzPreserved()
        {
            var path  = TempFile("simple.parquet");
            var entry = TestData.SimpleEntry();
            DiaNNParquetSpectralLibrary.Write(path, new[] { entry });

            var loaded = DiaNNParquetSpectralLibrary.Read(path);
            Assert.That(loaded.Count, Is.EqualTo(1));
            Assert.That(loaded[0].PrecursorMz, Is.EqualTo(entry.PrecursorMz).Within(1e-3));
        }

        [Test]
        public void WriteRead_SimpleEntry_FragmentCountPreserved()
        {
            var path  = TempFile("fragcount.parquet");
            var entry = TestData.SimpleEntry();
            DiaNNParquetSpectralLibrary.Write(path, new[] { entry });

            var loaded = DiaNNParquetSpectralLibrary.Read(path);
            Assert.That(loaded[0].Fragments.Count, Is.EqualTo(entry.Fragments.Count));
        }

        [Test]
        public void WriteRead_SmallLibrary_AllPrecursorsRoundTrip()
        {
            var path    = TempFile("library.parquet");
            var entries = TestData.SmallLibrary();
            DiaNNParquetSpectralLibrary.Write(path, entries);

            var loaded = DiaNNParquetSpectralLibrary.Read(path);
            Assert.That(loaded.Count, Is.EqualTo(entries.Count));
        }

        [Test]
        public void WriteRead_ModifiedSequences_Preserved()
        {
            var path  = TempFile("mods.parquet");
            var entry = TestData.OxidizedEntry();
            DiaNNParquetSpectralLibrary.Write(path, new[] { entry });

            var loaded = DiaNNParquetSpectralLibrary.Read(path);
            Assert.That(loaded[0].ModifiedSequence, Is.EqualTo(entry.ModifiedSequence));
        }

        [Test]
        public void WriteRead_RetentionTimesPreserved()
        {
            var path    = TempFile("rt.parquet");
            var entries = TestData.SmallLibrary().Take(4).ToList();
            DiaNNParquetSpectralLibrary.Write(path, entries);

            var loaded = DiaNNParquetSpectralLibrary.Read(path)
                .OrderBy(e => e.PrecursorMz).ToList();
            var original = entries.OrderBy(e => e.PrecursorMz).ToList();

            for (int i = 0; i < original.Count; i++)
                Assert.That(loaded[i].RetentionTime,
                    Is.EqualTo(original[i].RetentionTime).Within(1e-3),
                    $"RT mismatch at index {i}");
        }

        [Test]
        public void WriteRead_DecoyFlag_Preserved()
        {
            var path    = TempFile("decoy.parquet");
            var entries = new List<DiaNNLibraryEntry> { TestData.SimpleEntry(), TestData.DecoyEntry() };
            DiaNNParquetSpectralLibrary.Write(path, entries);

            var loaded = DiaNNParquetSpectralLibrary.Read(path);
            Assert.That(loaded.Any(e => e.IsDecoy),  Is.True);
            Assert.That(loaded.Any(e => !e.IsDecoy), Is.True);
        }

        [Test]
        public void WriteRead_IonMobility_Preserved()
        {
            var path  = TempFile("im.parquet");
            var entry = TestData.SimpleEntry();
            entry.IonMobility = 0.923;
            DiaNNParquetSpectralLibrary.Write(path, new[] { entry });

            var loaded = DiaNNParquetSpectralLibrary.Read(path);
            Assert.That(loaded[0].IonMobility, Is.EqualTo(entry.IonMobility).Within(1e-3));
        }

        [Test]
        public void WriteRead_NeutralLossTypes_Preserved()
        {
            var path  = TempFile("loss.parquet");
            var entry = TestData.OxidizedEntry(); // contains H2O loss
            DiaNNParquetSpectralLibrary.Write(path, new[] { entry });

            var loaded = DiaNNParquetSpectralLibrary.Read(path);
            Assert.That(loaded[0].Fragments.Any(f => f.LossType == "H2O"), Is.True);
        }

        [Test]
        public void ValidateSchema_OurOwnOutput_HasNoViolations()
        {
            var path = TempFile("validate.parquet");
            DiaNNParquetSpectralLibrary.Write(path, new[] { TestData.SimpleEntry() });

            var violations = DiaNNParquetSpectralLibrary.ValidateSchema(path);
            Assert.That(violations, Is.Empty,
                "Our own writer produced schema violations:\n" + string.Join("\n", violations));
        }

        [Test]
        public void GetSchema_FloatColumnsAreFloat_NotDouble()
        {
            var schema = DiaNNParquetSpectralLibrary.GetSchema();
            var floatCols = new[] { "PrecursorMz", "Tr_recalibrated", "IonMobility", "FragmentMz", "RelativeIntensity" };

            foreach (var colName in floatCols)
            {
                var field = schema.Fields.OfType<Parquet.Data.DataField>().FirstOrDefault(f => f.Name == colName);
                Assert.That(field, Is.Not.Null, $"Column '{colName}' not found in schema");
                Assert.That(field!.ClrType, Is.EqualTo(typeof(float)),
                    $"Column '{colName}' should be float but is {field.ClrType.Name}. DIA-NN requires FLOAT.");
            }
        }

        [Test]
        public void GetSchema_IntColumnsAreInt64_NotInt32()
        {
            var schema = DiaNNParquetSpectralLibrary.GetSchema();
            var int64Cols = new[] { "PrecursorCharge", "FragmentNumber", "FragmentCharge", "Proteotypic", "Decoy" };

            foreach (var colName in int64Cols)
            {
                var field = schema.Fields.OfType<Parquet.Data.DataField>().FirstOrDefault(f => f.Name == colName);
                Assert.That(field, Is.Not.Null, $"Column '{colName}' not found in schema");
                Assert.That(field!.ClrType, Is.EqualTo(typeof(long)),
                    $"Column '{colName}' should be long/INT64 but is {field.ClrType.Name}. DIA-NN requires INT64.");
            }
        }

        [Test]
        public void Write_EmptyList_ThrowsArgumentException()
        {
            var path = TempFile("empty.parquet");
            Assert.Throws<ArgumentException>(() =>
                DiaNNParquetSpectralLibrary.Write(path, new List<DiaNNLibraryEntry>()));
        }

        [Test]
        public void WriteInRowGroups_LargerLibrary_SameResultAsWrite()
        {
            var entries  = TestData.SmallLibrary();
            var pathA    = TempFile("single_group.parquet");
            var pathB    = TempFile("multi_group.parquet");

            DiaNNParquetSpectralLibrary.Write(pathA, entries);
            // Use rowsPerGroup = 3 to force multiple row groups on our 8-entry library
            DiaNNParquetSpectralLibrary.WriteInRowGroups(pathB, entries, rowsPerGroup: 3);

            var loadedA = DiaNNParquetSpectralLibrary.Read(pathA)
                .OrderBy(e => e.PrecursorMz).ToList();
            var loadedB = DiaNNParquetSpectralLibrary.Read(pathB)
                .OrderBy(e => e.PrecursorMz).ToList();

            Assert.That(loadedB.Count, Is.EqualTo(loadedA.Count));
            for (int i = 0; i < loadedA.Count; i++)
                Assert.That(loadedB[i].PrecursorMz, Is.EqualTo(loadedA[i].PrecursorMz).Within(1e-3));
        }

        [Test]
        public void WriteReadAsLibrarySpectra_PrecursorCountPreserved()
        {
            var path    = TempFile("spectra.parquet");
            var entries = TestData.SmallLibrary().Take(3).ToList();
            DiaNNParquetSpectralLibrary.Write(path, entries);

            var spectra = DiaNNParquetSpectralLibrary.ReadAsLibrarySpectra(path);
            Assert.That(spectra.Count, Is.EqualTo(entries.Count));
        }
    }

    // ═════════════════════════════════════════════════════════════════════════════════════════
    // 8. Integration tests: cross-format round-trips
    // ═════════════════════════════════════════════════════════════════════════════════════════

    [TestFixture]
    [Category("Integration")]
    public class DiaNNIntegrationTests
    {
        private string _tempDir = null!;

        [SetUp]
        public void SetUp() => _tempDir = Path.Combine(Path.GetTempPath(), "DiaNNIntegrationTests_" + Guid.NewGuid());
        [TearDown]
        public void TearDown() { if (Directory.Exists(_tempDir)) Directory.Delete(_tempDir, recursive: true); }

        private string TempFile(string name) { Directory.CreateDirectory(_tempDir); return Path.Combine(_tempDir, name); }

        [Test]
        public void TsvToParquetToTsv_SmallLibrary_PrecursorCountConsistent()
        {
            var entries  = TestData.SmallLibrary();
            var tsvPath  = TempFile("src.tsv");
            var pqPath   = TempFile("converted.parquet");
            var tsv2Path = TempFile("restored.tsv");

            // TSV → Parquet
            DiaNNTsvSpectralLibrary.WriteTsv(tsvPath, entries);
            var fromTsv = DiaNNTsvSpectralLibrary.ReadTsv(tsvPath);
            DiaNNParquetSpectralLibrary.Write(pqPath, fromTsv);

            // Parquet → TSV
            var fromPq = DiaNNParquetSpectralLibrary.Read(pqPath);
            DiaNNTsvSpectralLibrary.WriteTsv(tsv2Path, fromPq);
            var restored = DiaNNTsvSpectralLibrary.ReadTsv(tsv2Path);

            Assert.That(restored.Count, Is.EqualTo(entries.Count));
        }

        [Test]
        public void EntryToLibrarySpectrumToEntry_FragmentCountPreserved()
        {
            var entry   = TestData.SimpleEntry();
            var spectrum = entry.ToLibrarySpectrum();
            var restored = DiaNNLibraryEntry.FromLibrarySpectrum(spectrum);

            Assert.That(restored.Fragments.Count, Is.EqualTo(entry.Fragments.Count));
        }

        [Test]
        public void BinaryToIndexToQuery_FindsExpectedPrecursors()
        {
            var entries = TestData.SmallLibrary();
            var path    = TempFile("index_test.speclib");
            DiaNNSpecLibWriter.Write(path, entries);

            using var reader = new DiaNNSpecLibReader(path);
            var indexEntries = reader.ReadPrecursorIndex();

            using var index = new DiaNNSpecLibIndex(
                indexEntries,
                entryLoader: srcIdx => reader.ReadEntry(srcIdx));

            // Query around the first entry's m/z (≈ 401.5)
            using var scope = index.GetPrecursorsInWindow(400.0, 405.0);
            Assert.That(scope.Count, Is.GreaterThanOrEqualTo(1));
        }

        [Test]
        public void ParquetToIndex_DdaLookup_FindsEntryBySequenceAndCharge()
        {
            var entries = TestData.SmallLibrary().Take(3).ToList();
            var path    = TempFile("dda_lookup.parquet");
            DiaNNParquetSpectralLibrary.Write(path, entries);

            var loaded = DiaNNParquetSpectralLibrary.Read(path);
            using var index = DiaNNSpecLibIndex.BuildFromEntries(loaded);

            var first = entries[0];
            bool found = index.TryGetSpectrum(
                first.ModifiedSequence!, first.PrecursorCharge, out var found_entry);

            Assert.That(found, Is.True);
            Assert.That(found_entry, Is.Not.Null);
            Assert.That(found_entry!.PrecursorMz, Is.EqualTo(first.PrecursorMz).Within(1e-3));
        }

        [Test]
        public void FromKoinaPredictions_ToParquet_SchemaValid()
        {
            // Simulates the Koina → FragmentIntensityModel → LibrarySpectrum → Parquet pipeline
            var koinaSpectra = TestData.SmallLibrary()
                .Select(e => e.ToLibrarySpectrum())
                .ToList();

            var path = TempFile("koina_output.parquet");
            DiaNNParquetSpectralLibrary.WriteFromLibrarySpectra(path, koinaSpectra);

            var violations = DiaNNParquetSpectralLibrary.ValidateSchema(path);
            Assert.That(violations, Is.Empty,
                "Schema violations found in Koina-pipeline output:\n" + string.Join("\n", violations));

            var loaded = DiaNNParquetSpectralLibrary.Read(path);
            Assert.That(loaded.Count, Is.EqualTo(koinaSpectra.Count));
        }

        [Test]
        public void ModificationMapping_DiaNNRoundTripThroughTsv_SequencePreserved()
        {
            var entry = TestData.OxidizedEntry();
            var path  = TempFile("mod_roundtrip.tsv");

            DiaNNTsvSpectralLibrary.WriteTsv(path, new[] { entry });
            var loaded = DiaNNTsvSpectralLibrary.ReadTsv(path);

            Assert.That(loaded[0].ModifiedSequence, Is.EqualTo(entry.ModifiedSequence));
        }

        [Test]
        public void AllFormats_SameLibrary_YieldConsistentFragmentCounts()
        {
            var entries  = TestData.SmallLibrary().Take(3).ToList();
            var tsvPath  = TempFile("all.tsv");
            var pqPath   = TempFile("all.parquet");
            var binPath  = TempFile("all.speclib");

            DiaNNTsvSpectralLibrary.WriteTsv(tsvPath, entries);
            DiaNNParquetSpectralLibrary.Write(pqPath, entries);
            DiaNNSpecLibWriter.Write(binPath, entries);

            var fromTsv = DiaNNTsvSpectralLibrary.ReadTsv(tsvPath);
            var fromPq  = DiaNNParquetSpectralLibrary.Read(pqPath);

            using var reader  = new DiaNNSpecLibReader(binPath);
            var fromBin = reader.ReadAllEntries()
                .OrderBy(e => e.PrecursorMz).ToList();

            var sortedTsv    = fromTsv.OrderBy(e => e.PrecursorMz).ToList();
            var sortedPq     = fromPq.OrderBy(e => e.PrecursorMz).ToList();
            var sortedOrig   = entries.OrderBy(e => e.PrecursorMz).ToList();

            Assert.That(fromTsv.Count, Is.EqualTo(sortedOrig.Count), "TSV count mismatch");
            Assert.That(fromPq.Count,  Is.EqualTo(sortedOrig.Count), "Parquet count mismatch");
            Assert.That(fromBin.Count, Is.EqualTo(sortedOrig.Count), "Binary count mismatch");

            for (int i = 0; i < sortedOrig.Count; i++)
            {
                int origCount = sortedOrig[i].Fragments.Count;
                Assert.That(sortedTsv[i].Fragments.Count, Is.EqualTo(origCount),
                    $"TSV fragment count mismatch at entry {i}");
                Assert.That(sortedPq[i].Fragments.Count,  Is.EqualTo(origCount),
                    $"Parquet fragment count mismatch at entry {i}");
                Assert.That(fromBin[i].Fragments.Count,   Is.EqualTo(origCount),
                    $"Binary fragment count mismatch at entry {i}");
            }
        }
    }
}
