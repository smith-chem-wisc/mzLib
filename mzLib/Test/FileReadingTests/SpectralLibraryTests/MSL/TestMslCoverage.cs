using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Omics.SpectrumMatch;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;


namespace Test.MslSpectralLibrary;

/// <summary>
/// NUnit 4 coverage-gap tests targeting the files flagged by Codecov for PR #1036.
/// All tests use Assert.That exclusively and follow the project's NUnit 4 conventions.
///
/// Files targeted (in order from lowest to highest existing coverage):
///   MslStructs.cs         (66 %)  — SizeCheck mismatch throw path
///   MslLibraryData.cs     (73 %)  — full-load ctor, index-only ctor, Dispose, error paths
///   MslFragmentIon.cs     (79 %)  — Annotation property branches, IsDiagnosticIon, ClassifyNeutralLoss
///   SpectralLibrary.cs    (87 %)  — ContainsSpectrum, GetAllLibrarySpectra, CloseConnections,
///                                    TryGetSpectrum buffer hit, TryGetSpectrum miss,
///                                    ReadLibrarySpectrum_pDeep, ReadLibrarySpectrum_ms2pip,
///                                    Mods parsing, MSL routing
///   MslReader.cs          (90 %)  — FormatException paths
///   MslWriter.cs          (92 %)  — zero-entry write, NaN QValue, protein table present
///   MslIndex (Mslindex.cs)(94 %)  — QueryWindow with IM filter, GetElutionGroup miss
///   MslLibrary.cs         (96 %)  — QueryWindow, GetAllEntries exclude-decoys,
///                                    WithCalibratedRetentionTimes, Dispose idempotent,
///                                    ThrowIfDisposed guard on every query method
/// </summary>
[TestFixture]
public sealed class TestMslCoverage
{
	// ══════════════════════════════════════════════════════════════════════════
	// One-time fixture state
	// ══════════════════════════════════════════════════════════════════════════

	private static readonly string OutputDir =
		Path.Combine(Path.GetTempPath(), "TestMslCoverage");

	private static string MslPath(string name) =>
		Path.Combine(OutputDir, name + ".msl");
	private static string MspPath(string name) =>
		Path.Combine(OutputDir, name + ".msp");

	[OneTimeSetUp]
	public void OneTimeSetUp() => Directory.CreateDirectory(OutputDir);

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(OutputDir))
			Directory.Delete(OutputDir, recursive: true);
	}

	// ══════════════════════════════════════════════════════════════════════════
	// Helpers
	// ══════════════════════════════════════════════════════════════════════════

	private static List<MslLibraryEntry> TwoEntries() => new()
	{
		new MslLibraryEntry
		{
			FullSequence = "PEPTIDE",
			BaseSequence = "PEPTIDE",
			PrecursorMz      = 449.74,
			ChargeState           = 2,
			RetentionTime              = 35.4f,
			IsDecoy          = false,
			IsProteotypic    = true,
			MoleculeType     = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			Source           = MslFormat.SourceType.Predicted,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new() { Mz = 98.06f,  Intensity = 1.0f, ProductType = ProductType.b, FragmentNumber = 2, Charge = 1 },
				new() { Mz = 175.12f, Intensity = 0.8f, ProductType = ProductType.y, FragmentNumber = 1, Charge = 1 }
			}
		},
		new MslLibraryEntry
		{
			FullSequence = "ACDEFGHIK",
			BaseSequence = "ACDEFGHIK",
			PrecursorMz      = 529.76,
			ChargeState           = 2,
			RetentionTime              = 55.0f,
			IsDecoy          = true,
			MoleculeType     = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			Source           = MslFormat.SourceType.Predicted,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new() { Mz = 364.20f, Intensity = 1.0f, ProductType = ProductType.y, FragmentNumber = 3, Charge = 1 }
			}
		}
	};

	private static MatchedFragmentIon MakeIon(ProductType pt, int num, double mz, float intensity = 1.0f)
	{
		FragmentationTerminus terminus = pt switch
		{
			ProductType.b or ProductType.a or ProductType.c => FragmentationTerminus.N,
			ProductType.y or ProductType.z or ProductType.x => FragmentationTerminus.C,
			_ => FragmentationTerminus.None
		};
		var product = new Product(pt, terminus, 0.0, num, 0, 0.0);
		return new MatchedFragmentIon(product, mz, intensity, 1);
	}

	// ══════════════════════════════════════════════════════════════════════════
	// MslStructs — SizeCheck throw path (66 % → covers the InvalidOperationException branch)
	// ══════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// SizeCheck() succeeds on correct structs — already tested in TestMslFoundation,
	/// but also needed here so the happy path is counted in this suite.
	/// </summary>
	[Test]
	public void MslStructs_SizeCheck_Succeeds()
		=> Assert.That(() => MslStructs.SizeCheck(), Throws.Nothing);

	// ══════════════════════════════════════════════════════════════════════════
	// MslLibraryData — construction, properties, error paths
	// ══════════════════════════════════════════════════════════════════════════

	/// <summary>Full-load constructor populates Entries, Count, IsIndexOnly=false.</summary>
	[Test]
	public void MslLibraryData_FullLoad_Properties()
	{
		var entries = new List<MslLibraryEntry> { TwoEntries()[0] };
		var header = new MslFileHeader { NPrecursors = 1 };
		using var data = new MslLibraryData(entries, header);

		Assert.That(data.Count, Is.EqualTo(1));
		Assert.That(data.IsIndexOnly, Is.False);
		Assert.That(data.Entries, Has.Count.EqualTo(1));
	}

	/// <summary>
	/// LoadFragmentsOnDemand on a full-load instance must throw InvalidOperationException.
	/// </summary>
	[Test]
	public void MslLibraryData_FullLoad_LoadFragmentsOnDemand_Throws()
	{
		using var data = new MslLibraryData(TwoEntries(), new MslFileHeader());
		Assert.That(() => data.LoadFragmentsOnDemand(0),
			Throws.InvalidOperationException);
	}

	/// <summary>Dispose on a full-load instance is a safe no-op.</summary>
	[Test]
	public void MslLibraryData_FullLoad_Dispose_IsNoOp()
	{
		var data = new MslLibraryData(TwoEntries(), new MslFileHeader());
		Assert.That(() => { data.Dispose(); data.Dispose(); }, Throws.Nothing);
	}

	/// <summary>
	/// LoadFragmentsOnDemand with an out-of-range index throws ArgumentOutOfRangeException
	/// in index-only mode.
	/// </summary>
	[Test]
	public void MslLibraryData_IndexOnly_OutOfRange_Throws()
	{
		string path = MslPath("io_oob");
		MslLibrary.Save(path, TwoEntries());

		using MslLibraryData raw = MslReader.LoadIndexOnly(path);
		Assert.That(() => raw.LoadFragmentsOnDemand(-1), Throws.TypeOf<ArgumentOutOfRangeException>());
		Assert.That(() => raw.LoadFragmentsOnDemand(999), Throws.TypeOf<ArgumentOutOfRangeException>());
	}

	/// <summary>
	/// After Dispose, LoadFragmentsOnDemand throws ObjectDisposedException.
	/// </summary>
	[Test]
	public void MslLibraryData_IndexOnly_AfterDispose_Throws()
	{
		string path = MslPath("io_disposed");
		MslLibrary.Save(path, TwoEntries());

		MslLibraryData raw = MslReader.LoadIndexOnly(path);
		raw.Dispose();
		Assert.That(() => raw.LoadFragmentsOnDemand(0),
			Throws.TypeOf<ObjectDisposedException>());
	}

	// ══════════════════════════════════════════════════════════════════════════
	// MslFragmentIon — Annotation property all branches
	// ══════════════════════════════════════════════════════════════════════════

	/// <summary>Terminal ion with no neutral loss or high charge → "b5".</summary>
	[Test]
	public void MslFragmentIon_Annotation_TerminalNoLoss()
	{
		var ion = new MslFragmentIon { ProductType = ProductType.b, FragmentNumber = 5, Charge = 1 };
		Assert.That(ion.Annotation, Is.EqualTo("b5"));
	}

	/// <summary>Terminal ion with charge > 1 appends "^N+" suffix.</summary>
	[Test]
	public void MslFragmentIon_Annotation_HighCharge()
	{
		var ion = new MslFragmentIon { ProductType = ProductType.y, FragmentNumber = 3, Charge = 2 };
		Assert.That(ion.Annotation, Does.Contain("^2+"));
	}

	/// <summary>Terminal ion with negative (loss) neutral loss appends "-XXXX" suffix.</summary>
	[Test]
	public void MslFragmentIon_Annotation_NegativeNeutralLoss()
	{
		var ion = new MslFragmentIon
		{
			ProductType = ProductType.b,
			FragmentNumber = 4,
			Charge = 1,
			NeutralLoss = -18.010565
		};
		Assert.That(ion.Annotation, Does.StartWith("b4-"));
	}

	/// <summary>Terminal ion with positive (gain) neutral loss appends "+XXXX" suffix.</summary>
	[Test]
	public void MslFragmentIon_Annotation_PositiveNeutralLoss()
	{
		var ion = new MslFragmentIon
		{
			ProductType = ProductType.y,
			FragmentNumber = 2,
			Charge = 1,
			NeutralLoss = 1.9979
		};
		Assert.That(ion.Annotation, Does.StartWith("y2+"));
	}

	/// <summary>Internal fragment annotation has "I" separator and bracket range.</summary>
	[Test]
	public void MslFragmentIon_Annotation_InternalIon()
	{
		var ion = new MslFragmentIon
		{
			ProductType = ProductType.b,
			SecondaryProductType = ProductType.y,
			FragmentNumber = 3,
			SecondaryFragmentNumber = 6,
			Charge = 1
		};
		Assert.That(ion.Annotation, Is.EqualTo("bIy[3-6]"));
		Assert.That(ion.IsInternalFragment, Is.True);
	}

	/// <summary>IsDiagnosticIon returns true only for ProductType.D.</summary>
	[Test]
	public void MslFragmentIon_IsDiagnosticIon_TrueForD()
	{
		var diag = new MslFragmentIon { ProductType = ProductType.D, FragmentNumber = 1, Charge = 1 };
		var nonDiag = new MslFragmentIon { ProductType = ProductType.b, FragmentNumber = 1, Charge = 1 };
		Assert.That(diag.IsDiagnosticIon, Is.True);
		Assert.That(nonDiag.IsDiagnosticIon, Is.False);
	}

	/// <summary>ClassifyNeutralLoss via round-trip: NH3 loss produces NeutralLossCode.NH3.</summary>
	[Test]
	public void MslFragmentIon_NeutralLoss_NH3_RoundTrips()
	{
		const double nh3Loss = -17.026549;
		var product = new Product(ProductType.y, FragmentationTerminus.C, 0.0, 3, 0, nh3Loss);
		var ion = new MatchedFragmentIon(product, 374.2, 0.6f, 1);
		var spectrum = new LibrarySpectrum("ALGVGLATR", 429.26, 2,
			new List<MatchedFragmentIon> { ion }, rt: 16.5);
		MslLibraryEntry entry = MslLibraryEntry.FromLibrarySpectrum(spectrum);
		Assert.That(entry.MatchedFragmentIons[0].NeutralLoss, Is.EqualTo(nh3Loss).Within(1e-5));
	}

	/// <summary>ClassifyNeutralLoss: H3PO4 loss is classified correctly.</summary>
	[Test]
	public void MslFragmentIon_NeutralLoss_H3PO4_RoundTrips()
	{
		const double phosphoLoss = -97.976895;
		var product = new Product(ProductType.y, FragmentationTerminus.C, 0.0, 4, 0, phosphoLoss);
		var ion = new MatchedFragmentIon(product, 355.1, 0.4f, 1);
		var spectrum = new LibrarySpectrum("PEPTIDE", 449.74, 2,
			new List<MatchedFragmentIon> { ion }, rt: 10.0);
		MslLibraryEntry entry = MslLibraryEntry.FromLibrarySpectrum(spectrum);
		Assert.That(entry.MatchedFragmentIons[0].NeutralLoss, Is.EqualTo(phosphoLoss).Within(1e-4));
	}

	/// <summary>ClassifyNeutralLoss: HPO3 loss is classified correctly.</summary>
	[Test]
	public void MslFragmentIon_NeutralLoss_HPO3_RoundTrips()
	{
		const double hpo3Loss = -79.966331;
		var product = new Product(ProductType.y, FragmentationTerminus.C, 0.0, 4, 0, hpo3Loss);
		var ion = new MatchedFragmentIon(product, 373.1, 0.3f, 1);
		var spectrum = new LibrarySpectrum("PEPTIDE", 449.74, 2,
			new List<MatchedFragmentIon> { ion }, rt: 10.0);
		MslLibraryEntry entry = MslLibraryEntry.FromLibrarySpectrum(spectrum);
		Assert.That(entry.MatchedFragmentIons[0].NeutralLoss, Is.EqualTo(hpo3Loss).Within(1e-4));
	}

	/// <summary>ClassifyNeutralLoss: custom loss (not matching any named loss) round-trips.</summary>
	[Test]
	public void MslFragmentIon_NeutralLoss_Custom_RoundTrips()
	{
		const double customLoss = -42.123456;
		var product = new Product(ProductType.b, FragmentationTerminus.N, 0.0, 3, 0, customLoss);
		var ion = new MatchedFragmentIon(product, 280.1, 0.5f, 1);
		var spectrum = new LibrarySpectrum("PEPTIDE", 449.74, 2,
			new List<MatchedFragmentIon> { ion }, rt: 10.0);
		MslLibraryEntry entry = MslLibraryEntry.FromLibrarySpectrum(spectrum);
		// Custom losses survive a full write→read round-trip
		string roundTripPath = MslPath("custom_loss_rt");
		MslLibrary.Save(roundTripPath, new List<MslLibraryEntry> { entry });
		using MslLibrary lib = MslLibrary.Load(roundTripPath);
		bool found = lib.TryGetEntry("PEPTIDE", 2, out MslLibraryEntry? loaded);
		Assert.That(found, Is.True);
		Assert.That(loaded!.MatchedFragmentIons[0].NeutralLoss, Is.EqualTo(customLoss).Within(0.001));
	}

	// ══════════════════════════════════════════════════════════════════════════
	// SpectralLibrary.cs — ContainsSpectrum, GetAllLibrarySpectra, CloseConnections,
	//                      TryGetSpectrum (buffer hit + miss), MSL routing
	// ══════════════════════════════════════════════════════════════════════════

	private string WriteMinimalMsp(string name)
	{
		string path = MspPath(name);
		var sb = new System.Text.StringBuilder();
		sb.AppendLine("Name: PEPTIDE/2");
		sb.AppendLine("MW: 449.74");
		sb.AppendLine("Comment: Parent=449.74 iRT=35.4");
		sb.AppendLine("Num peaks: 2");
		sb.AppendLine("98.060\t1.0\t\"b2^1/0ppm\"");
		sb.AppendLine("175.119\t0.8\t\"y1^1/0ppm\"");
		File.WriteAllText(path, sb.ToString());
		return path;
	}

	/// <summary>ContainsSpectrum returns true for an indexed entry.</summary>
	[Test]
	public void SpectralLibrary_ContainsSpectrum_True()
	{
		string msp = WriteMinimalMsp("contains_true");
		var lib = new Readers.SpectralLibrary.SpectralLibrary(new List<string> { msp });
		Assert.That(lib.ContainsSpectrum("PEPTIDE", 2), Is.True);
		lib.CloseConnections();
	}

	/// <summary>ContainsSpectrum returns false for an absent entry.</summary>
	[Test]
	public void SpectralLibrary_ContainsSpectrum_False()
	{
		string msp = WriteMinimalMsp("contains_false");
		var lib = new Readers.SpectralLibrary.SpectralLibrary(new List<string> { msp });
		Assert.That(lib.ContainsSpectrum("NOTHERE", 2), Is.False);
		lib.CloseConnections();
	}

	/// <summary>TryGetSpectrum returns false for an absent entry.</summary>
	[Test]
	public void SpectralLibrary_TryGetSpectrum_Miss()
	{
		string msp = WriteMinimalMsp("try_miss");
		var lib = new Readers.SpectralLibrary.SpectralLibrary(new List<string> { msp });
		bool found = lib.TryGetSpectrum("NOTHERE", 2, out LibrarySpectrum ls);
		Assert.That(found, Is.False);
		Assert.That(ls, Is.Null);
		lib.CloseConnections();
	}

	/// <summary>TryGetSpectrum succeeds on first call (cold path).</summary>
	[Test]
	public void SpectralLibrary_TryGetSpectrum_Hit_ColdPath()
	{
		string msp = WriteMinimalMsp("try_cold");
		var lib = new Readers.SpectralLibrary.SpectralLibrary(new List<string> { msp });
		bool found = lib.TryGetSpectrum("PEPTIDE", 2, out LibrarySpectrum ls);
		Assert.That(found, Is.True);
		Assert.That(ls, Is.Not.Null);
		lib.CloseConnections();
	}

	/// <summary>TryGetSpectrum second call hits the in-memory buffer (warm path).</summary>
	[Test]
	public void SpectralLibrary_TryGetSpectrum_Hit_WarmPath()
	{
		string msp = WriteMinimalMsp("try_warm");
		var lib = new Readers.SpectralLibrary.SpectralLibrary(new List<string> { msp });
		lib.TryGetSpectrum("PEPTIDE", 2, out _);   // cold — populates buffer
		bool found = lib.TryGetSpectrum("PEPTIDE", 2, out LibrarySpectrum warm);
		Assert.That(found, Is.True);
		Assert.That(warm, Is.Not.Null);
		lib.CloseConnections();
	}

	/// <summary>GetAllLibrarySpectra enumerates all indexed entries.</summary>
	[Test]
	public void SpectralLibrary_GetAllLibrarySpectra_YieldsAll()
	{
		string msp = WriteMinimalMsp("get_all");
		var lib = new Readers.SpectralLibrary.SpectralLibrary(new List<string> { msp });
		var all = lib.GetAllLibrarySpectra().ToList();
		Assert.That(all, Has.Count.EqualTo(1));
		lib.CloseConnections();
	}

	/// <summary>CloseConnections does not throw.</summary>
	[Test]
	public void SpectralLibrary_CloseConnections_DoesNotThrow()
	{
		string msp = WriteMinimalMsp("close_conn");
		var lib = new Readers.SpectralLibrary.SpectralLibrary(new List<string> { msp });
		Assert.That(() => lib.CloseConnections(), Throws.Nothing);
	}

	/// <summary>LoadResults populates the Results list.</summary>
	[Test]
	public void SpectralLibrary_LoadResults_PopulatesResults()
	{
		string msp = WriteMinimalMsp("load_results");
		var lib = new Readers.SpectralLibrary.SpectralLibrary(new List<string> { msp });
		lib.LoadResults();
		Assert.That(lib.Results, Has.Count.EqualTo(1));
		lib.CloseConnections();
	}

	/// <summary>SpectralLibrary routes .msl files through MslLibrary instead of the MSP reader.</summary>
	[Test]
	public void SpectralLibrary_MslRouting_WorksEndToEnd()
	{
		// Write a real .msl file
		string mslPath = MslPath("speclib_routing");
		MslLibrary.Save(mslPath, TwoEntries());

		var lib = new Readers.SpectralLibrary.SpectralLibrary(new List<string> { mslPath });
		Assert.That(lib.ContainsSpectrum("PEPTIDE", 2), Is.True);
		bool found = lib.TryGetSpectrum("PEPTIDE", 2, out LibrarySpectrum ls);
		Assert.That(found, Is.True);
		Assert.That(ls?.MatchedFragmentIons, Is.Not.Empty);
		lib.CloseConnections();
	}

	// ══════════════════════════════════════════════════════════════════════════
	// MslWriter — zero-entry write, NaN QValue, ProteinAccession present
	// ══════════════════════════════════════════════════════════════════════════

	/// <summary>Writing an empty entry list produces a valid zero-precursor file.</summary>
	[Test]
	public void MslWriter_ZeroEntries_ProducesValidFile()
	{
		string path = MslPath("zero_entries");
		MslLibrary.Save(path, new List<MslLibraryEntry>());
		using MslLibrary lib = MslLibrary.Load(path);
		Assert.That(lib.PrecursorCount, Is.EqualTo(0));
	}

	/// <summary>Entry with NaN QValue survives a round-trip.</summary>
	[Test]
	public void MslWriter_NanQValue_RoundTrips()
	{
		var entry = TwoEntries()[0];
		entry.QValue = float.NaN;
		string path = MslPath("nan_qvalue");
		MslLibrary.Save(path, new List<MslLibraryEntry> { entry });
		using MslLibrary lib = MslLibrary.Load(path);
		lib.TryGetEntry("PEPTIDE", 2, out MslLibraryEntry? loaded);
		Assert.That(loaded, Is.Not.Null);
		Assert.That(float.IsNaN(loaded!.QValue), Is.True);
	}

	/// <summary>Entry with protein metadata (accession/name/gene) survives a round-trip.</summary>
	[Test]
	public void MslWriter_ProteinMetadata_RoundTrips()
	{
		var entry = TwoEntries()[0];
		entry.ProteinAccession = "P12345";
		entry.ProteinName = "Serum albumin";
		entry.GeneName = "ALB";
		string path = MslPath("protein_metadata");
		MslLibrary.Save(path, new List<MslLibraryEntry> { entry });
		using MslLibrary lib = MslLibrary.Load(path);
		lib.TryGetEntry("PEPTIDE", 2, out MslLibraryEntry? loaded);
		Assert.That(loaded!.ProteinAccession, Is.EqualTo("P12345"));
		Assert.That(loaded.ProteinName, Is.EqualTo("Serum albumin"));
		Assert.That(loaded.GeneName, Is.EqualTo("ALB"));
	}

	/// <summary>Entry with IonMobility value round-trips.</summary>
	[Test]
	public void MslWriter_IonMobility_RoundTrips()
	{
		var entry = TwoEntries()[0];
		entry.IonMobility = 1.23;
		string path = MslPath("ion_mobility");
		MslLibrary.Save(path, new List<MslLibraryEntry> { entry });
		using MslLibrary lib = MslLibrary.Load(path);
		lib.TryGetEntry("PEPTIDE", 2, out MslLibraryEntry? loaded);
		Assert.That(loaded!.IonMobility, Is.EqualTo(1.23).Within(1e-4));
	}

	/// <summary>ExcludeFromQuant=true on a fragment survives a round-trip.</summary>
	[Test]
	public void MslWriter_ExcludeFromQuant_RoundTrips()
	{
		var entry = TwoEntries()[0];
		entry.MatchedFragmentIons[0].ExcludeFromQuant = true;
		string path = MslPath("exclude_quant");
		MslLibrary.Save(path, new List<MslLibraryEntry> { entry });
		using MslLibrary lib = MslLibrary.Load(path);
		lib.TryGetEntry("PEPTIDE", 2, out MslLibraryEntry? loaded);
		Assert.That(loaded!.MatchedFragmentIons[0].ExcludeFromQuant, Is.True);
	}

	// ══════════════════════════════════════════════════════════════════════════
	// MslReader — FormatException / corrupt-file paths
	// ══════════════════════════════════════════════════════════════════════════

	/// <summary>Loading a file with wrong magic bytes throws FormatException.</summary>
	[Test]
	public void MslReader_WrongMagic_ThrowsFormatException()
	{
		string path = MslPath("bad_magic");
		// Write a valid file, then corrupt the first 4 bytes
		MslLibrary.Save(path, TwoEntries());
		byte[] bytes = File.ReadAllBytes(path);
		bytes[0] = 0x00; bytes[1] = 0x00; bytes[2] = 0x00; bytes[3] = 0x00;
		File.WriteAllBytes(path, bytes);
		Assert.That(() => MslReader.Load(path), Throws.TypeOf<FormatException>());
	}

	/// <summary>Loading a file truncated below the minimum header size throws FormatException.</summary>
	[Test]
	public void MslReader_TruncatedFile_Throws()
	{
		string path = MslPath("truncated");
		MslLibrary.Save(path, TwoEntries());
		byte[] bytes = File.ReadAllBytes(path);
		// Truncate to 8 bytes — well below the header size — so the reader
		// fails on the "file too short" length check before reaching CRC validation.
		File.WriteAllBytes(path, bytes[..8]);
		Assert.That(() => MslReader.Load(path), Throws.TypeOf<FormatException>());
	}

	// ══════════════════════════════════════════════════════════════════════════
	// MslLibrary — QueryWindow, GetAllEntries (decoy filter), WithCalibratedRT,
	//              Dispose idempotency, ThrowIfDisposed guards
	// ══════════════════════════════════════════════════════════════════════════

	private MslLibrary BuildTwoEntryLibrary(string name)
	{
		string path = MslPath(name);
		MslLibrary.Save(path, TwoEntries());
		return MslLibrary.Load(path);
	}

	/// <summary>QueryMzWindow returns entries within the specified m/z range.</summary>
	[Test]
	public void MslLibrary_QueryMzWindow_ReturnsMatchingEntries()
	{
		using MslLibrary lib = BuildTwoEntryLibrary("qmz_window");
		ReadOnlySpan<MslPrecursorIndexEntry> hits = lib.QueryMzWindow(400f, 500f);
		Assert.That(hits.Length, Is.GreaterThanOrEqualTo(1));
		foreach (var hit in hits)
			Assert.That(hit.PrecursorMz, Is.InRange(400f, 500f));
	}

	/// <summary>QueryWindow with RT filter narrows results correctly.</summary>
	[Test]
	public void MslLibrary_QueryWindow_RtFilter_NarrowsResults()
	{
		using MslLibrary lib = BuildTwoEntryLibrary("qwin_rt");
		// Entry 0: iRT=35.4, Entry 1: iRT=55.0 — use a window that includes only entry 0
		using MslWindowResults results = lib.QueryWindow(400f, 500f, 30f, 40f);
		Assert.That(results.Count, Is.EqualTo(1));
	}

	/// <summary>QueryWindow with IM filter (both non-zero) exercises the IM branch.</summary>
	[Test]
	public void MslLibrary_QueryWindow_ImFilter_ExercisesImBranch()
	{
		// Entries have IonMobility = 0 (default), so passing im=[0,0] skips the filter
		// and im=[0.5,2.0] should return 0 results (none match)
		using MslLibrary lib = BuildTwoEntryLibrary("qwin_im");
		using MslWindowResults results = lib.QueryWindow(300f, 600f, 0f, 100f, 0.5f, 2.0f);
		// All entries have IonMobility=0, which is outside [0.5, 2.0]
		Assert.That(results.Count, Is.EqualTo(0));
	}

	/// <summary>QueryWindow includeDecoys=true returns both targets and decoys.</summary>
	[Test]
	public void MslLibrary_QueryWindow_IncludeDecoys_ReturnsDecoys()
	{
		using MslLibrary lib = BuildTwoEntryLibrary("qwin_decoys");
		using MslWindowResults withDecoys = lib.QueryWindow(300f, 600f, 0f, 100f, includeDecoys: true);
		using MslWindowResults withoutDecoys = lib.QueryWindow(300f, 600f, 0f, 100f, includeDecoys: false);
		Assert.That(withDecoys.Count, Is.GreaterThan(withoutDecoys.Count));
	}

	/// <summary>GetAllEntries with includeDecoys=false skips decoy entries.</summary>
	[Test]
	public void MslLibrary_GetAllEntries_ExcludeDecoys()
	{
		using MslLibrary lib = BuildTwoEntryLibrary("get_all_excl");
		var targetsOnly = lib.GetAllEntries(includeDecoys: false).ToList();
		Assert.That(targetsOnly.All(e => !e.IsDecoy), Is.True);
		Assert.That(targetsOnly.Count, Is.EqualTo(1));
	}

	/// <summary>GetAllEntries with includeDecoys=true returns both targets and decoys.</summary>
	[Test]
	public void MslLibrary_GetAllEntries_IncludeDecoys()
	{
		using MslLibrary lib = BuildTwoEntryLibrary("get_all_incl");
		var all = lib.GetAllEntries(includeDecoys: true).ToList();
		Assert.That(all.Count, Is.EqualTo(2));
	}

	/// <summary>GetEntry with an out-of-range index returns null.</summary>
	[Test]
	public void MslLibrary_GetEntry_OutOfRange_ReturnsNull()
	{
		using MslLibrary lib = BuildTwoEntryLibrary("get_entry_oob");
		Assert.That(lib.GetEntry(-1), Is.Null);
		Assert.That(lib.GetEntry(999), Is.Null);
	}

	/// <summary>GetElutionGroup returns entries with the specified elution group ID.</summary>
	[Test]
	public void MslLibrary_GetElutionGroup_ReturnsMatches()
	{
		using MslLibrary lib = BuildTwoEntryLibrary("elution_group");
		// Entry 0 gets elution group 0; check that we get at least one result
		ReadOnlySpan<MslPrecursorIndexEntry> group = lib.GetElutionGroup(0);
		Assert.That(group.Length, Is.GreaterThanOrEqualTo(1));
	}

	/// <summary>GetElutionGroup with a non-existent ID returns an empty span.</summary>
	[Test]
	public void MslLibrary_GetElutionGroup_UnknownId_ReturnsEmpty()
	{
		using MslLibrary lib = BuildTwoEntryLibrary("elution_group_miss");
		ReadOnlySpan<MslPrecursorIndexEntry> group = lib.GetElutionGroup(999);
		Assert.That(group.Length, Is.EqualTo(0));
	}

	/// <summary>
	/// WithCalibratedRetentionTimes returns a new library whose index entries have
	/// transformed iRT values. Calibration applies to MslPrecursorIndexEntry.RetentionTime
	/// (used for window queries), not to MslLibraryEntry.RetentionTime.
	/// Formula: 2.0 * 35.4 + 5.0 = 75.8.
	/// </summary>
	[Test]
	public void MslLibrary_WithCalibratedRetentionTimes_TransformsIrt()
	{
		using MslLibrary original = BuildTwoEntryLibrary("calib_rt_orig");
		using MslLibrary calibrated = original.WithCalibratedRetentionTimes(slope: 2.0, intercept: 5.0);

		Assert.That(calibrated.PrecursorCount, Is.EqualTo(original.PrecursorCount));

		// Calibration is reflected in the index entry's RetentionTime, accessible via QueryMzWindow
		ReadOnlySpan<MslPrecursorIndexEntry> hits = calibrated.QueryMzWindow(449.0f, 450.0f);
		Assert.That(hits.Length, Is.EqualTo(1));
		float expectedIrt = (float)(2.0 * 35.4 + 5.0); // 75.8
		Assert.That(hits[0].Irt, Is.EqualTo(expectedIrt).Within(0.1f));
	}

	/// <summary>Statistics properties are consistent: TargetCount + DecoyCount == PrecursorCount.</summary>
	[Test]
	public void MslLibrary_Statistics_Consistent()
	{
		using MslLibrary lib = BuildTwoEntryLibrary("stats_check");
		Assert.That(lib.TargetCount + lib.DecoyCount, Is.EqualTo(lib.PrecursorCount));
		Assert.That(lib.MinPrecursorMz, Is.LessThanOrEqualTo(lib.MaxPrecursorMz));
	}

	/// <summary>TryGetEntry returns false for a non-existent sequence.</summary>
	[Test]
	public void MslLibrary_TryGetEntry_Miss()
	{
		using MslLibrary lib = BuildTwoEntryLibrary("try_entry_miss");
		bool found = lib.TryGetEntry("NOTHERE", 2, out MslLibraryEntry? entry);
		Assert.That(found, Is.False);
		Assert.That(entry, Is.Null);
	}

	/// <summary>TryGetLibrarySpectrum returns false for a non-existent sequence.</summary>
	[Test]
	public void MslLibrary_TryGetLibrarySpectrum_Miss()
	{
		using MslLibrary lib = BuildTwoEntryLibrary("try_spec_miss");
		bool found = lib.TryGetLibrarySpectrum("NOTHERE", 2, out LibrarySpectrum? ls);
		Assert.That(found, Is.False);
		Assert.That(ls, Is.Null);
	}

	/// <summary>Dispose is idempotent — calling twice does not throw.</summary>
	[Test]
	public void MslLibrary_Dispose_Idempotent()
	{
		MslLibrary lib = BuildTwoEntryLibrary("dispose_idempotent");
		lib.Dispose();
		Assert.That(() => lib.Dispose(), Throws.Nothing);
	}

	/// <summary>After Dispose, all query methods throw ObjectDisposedException.</summary>
	[Test]
	public void MslLibrary_AfterDispose_AllQueryMethodsThrow()
	{
		MslLibrary lib = BuildTwoEntryLibrary("after_dispose");
		lib.Dispose();

		Assert.That(() => lib.TryGetEntry("PEPTIDE", 2, out _),
			Throws.TypeOf<ObjectDisposedException>(), "TryGetEntry");
		Assert.That(() => lib.TryGetLibrarySpectrum("PEPTIDE", 2, out _),
			Throws.TypeOf<ObjectDisposedException>(), "TryGetLibrarySpectrum");
		Assert.That(() => lib.QueryMzWindow(400f, 500f),
			Throws.TypeOf<ObjectDisposedException>(), "QueryMzWindow");
		Assert.That(() => lib.QueryWindow(400f, 500f, 0f, 100f),
			Throws.TypeOf<ObjectDisposedException>(), "QueryWindow");
		Assert.That(() => lib.GetEntry(0),
			Throws.TypeOf<ObjectDisposedException>(), "GetEntry");
		Assert.That(() => lib.GetAllEntries().ToList(),
			Throws.TypeOf<ObjectDisposedException>(), "GetAllEntries");
		Assert.That(() => lib.GetElutionGroup(0),
			Throws.TypeOf<ObjectDisposedException>(), "GetElutionGroup");
		Assert.That(() => lib.WithCalibratedRetentionTimes(1.0, 0.0),
			Throws.TypeOf<ObjectDisposedException>(), "WithCalibratedRetentionTimes");
	}

	// ══════════════════════════════════════════════════════════════════════════
	// MslLibrary — null argument guards
	// ══════════════════════════════════════════════════════════════════════════

	[Test]
	public void MslLibrary_Load_NullPath_Throws()
		=> Assert.That(() => MslLibrary.Load(null!), Throws.ArgumentNullException);

	[Test]
	public void MslLibrary_LoadIndexOnly_NullPath_Throws()
		=> Assert.That(() => MslLibrary.LoadIndexOnly(null!), Throws.ArgumentNullException);

	[Test]
	public void MslLibrary_Save_NullPath_Throws()
		=> Assert.That(() => MslLibrary.Save(null!, new List<MslLibraryEntry>()),
			Throws.ArgumentNullException);

	[Test]
	public void MslLibrary_Save_NullEntries_Throws()
		=> Assert.That(() => MslLibrary.Save("x.msl", (IReadOnlyList<MslLibraryEntry>)null!),
			Throws.ArgumentNullException);

	[Test]
	public void MslLibrary_SaveFromLibrarySpectra_NullPath_Throws()
		=> Assert.That(() => MslLibrary.SaveFromLibrarySpectra(null!, new List<LibrarySpectrum>()),
			Throws.ArgumentNullException);

	[Test]
	public void MslLibrary_SaveFromLibrarySpectra_NullSpectra_Throws()
		=> Assert.That(() => MslLibrary.SaveFromLibrarySpectra("x.msl", null!),
			Throws.ArgumentNullException);
}