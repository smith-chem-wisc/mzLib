using Chemistry;
using Chromatography.RetentionTimeCalibration;
using MzLibUtil;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Omics.SpectrumMatch;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

namespace Test.MslSpectralLibrary;

/// <summary>
/// Tests for the .msl correctness fixes: the neutral-loss sign, the merger's choice of winner,
/// the window result's Dispose, unconvertible spectra, values the writer stored wrongly,
/// compression levels, temp-file names, disposal during enumeration, the calibration slope,
/// the header version check and the proton constant. Each fails on the code before its fix,
/// except one guard that pins behaviour the merger rewrite keeps: KeepLast keeps the last
/// file's copy.
/// </summary>
[TestFixture]
public class TestMslCorrectnessFixes
{
	private static string _dir = null!;

	[OneTimeSetUp]
	public void OneTimeSetUp()
	{
		_dir = Path.Combine(TestContext.CurrentContext.TestDirectory, nameof(TestMslCorrectnessFixes));
		Directory.CreateDirectory(_dir);
	}

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(_dir))
			Directory.Delete(_dir, recursive: true);
	}

	private static string PathFor(string name) => Path.Combine(_dir, name + ".msl");

	private static MslLibraryEntry Entry(string sequence, double mz, float rt = 10f, float qValue = float.NaN,
		double neutralLoss = 0.0, int fragmentCharge = 1)
	{
		return new MslLibraryEntry
		{
			FullSequence = sequence,
			BaseSequence = sequence,
			PrecursorMz = mz,
			ChargeState = 2,
			RetentionTime = rt,
			QValue = qValue,
			Source = MslFormat.SourceType.Empirical,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz = 500f, Intensity = 1f, ProductType = ProductType.y, FragmentNumber = 4,
					ResiduePosition = 4, Charge = fragmentCharge, NeutralLoss = neutralLoss
				}
			}
		};
	}

	private static LibrarySpectrum Spectrum(string sequence, int charge, double neutralLoss = 0.0,
		double? rt = 12.0, int fragments = 1, double precursorMz = 450.0)
	{
		var ions = new List<MatchedFragmentIon>();
		for (int i = 0; i < fragments; i++)
		{
			var product = new Product(ProductType.y, FragmentationTerminus.C, 600.0 + i, 5 + i, 5 + i, neutralLoss);
			ions.Add(new MatchedFragmentIon(product, 502.0 + i, 1.0, 1));
		}
		return new LibrarySpectrum(sequence, precursorMz, charge, ions, rt);
	}

	// ── Neutral-loss sign ────────────────────────────────────────────────────

	/// <summary>
	/// A loss as mzLib's Product carries it (positive, as the MSP reader and Koina produce it)
	/// is stored as its named code and reads back positive. Before the fix it was stored as a
	/// Custom value and annotated as a gain ("y5+97.9769").
	/// </summary>
	[Test]
	public void ProductLoss_IsStoredAsItsNamedCode_AndReadsBackAsALoss()
	{
		string path = PathFor(nameof(ProductLoss_IsStoredAsItsNamedCode_AndReadsBackAsALoss));
		MslWriter.WriteFromLibrarySpectra(path, new[] { Spectrum("PEPTSIDE", 2, neutralLoss: 97.976895) });

		Assert.That(MslReader.ReadHeaderOnly(path).FileFlags & MslFormat.FileFlagHasExtAnnotations, Is.EqualTo(0),
			"A named loss needs no custom-loss table.");
		using MslLibrary library = MslLibrary.Load(path);
		MslFragmentIon read = library.GetEntry(0, includeFragments: true)!.MatchedFragmentIons.Single();
		Assert.That(read.NeutralLoss, Is.EqualTo(97.976895).Within(1e-6));
		Assert.That(read.Annotation, Is.EqualTo("y5-97.9769"));
	}

	/// <summary>
	/// A named loss reaches a LibrarySpectrum (what MetaMorpheus reads) with Product's sign, so
	/// Product's annotation shows a loss. Before the fix the code decoded to -18.01, and the
	/// annotation read "y4--18.01".
	/// </summary>
	[Test]
	public void NamedLoss_ReachesLibrarySpectrum_AsAProductLoss()
	{
		string path = PathFor(nameof(NamedLoss_ReachesLibrarySpectrum_AsAProductLoss));
		MslWriter.Write(path, new[] { Entry("PEPTIDEK", 450.0, neutralLoss: 18.010565) });

		Assert.That(MslReader.ReadHeaderOnly(path).FileFlags & MslFormat.FileFlagHasExtAnnotations, Is.EqualTo(0),
			"Water is a named loss.");
		using MslLibrary library = MslLibrary.Load(path);
		Assert.That(library.TryGetLibrarySpectrum("PEPTIDEK", 2, out LibrarySpectrum? spectrum), Is.True);
		Product product = spectrum!.MatchedFragmentIons.Single().NeutralTheoreticalProduct;
		Assert.That(product.NeutralLoss, Is.EqualTo(18.010565).Within(1e-6));
		Assert.That(product.Annotation, Is.EqualTo("y4-18.01"));
	}

	/// <summary>
	/// A negative value is a gain: it is kept exactly as a Custom value, sign and all. Before
	/// the fix -18.01 was taken for the named water loss.
	/// </summary>
	[Test]
	public void NegativeLoss_IsAGain_KeptAsACustomValue()
	{
		string path = PathFor(nameof(NegativeLoss_IsAGain_KeptAsACustomValue));
		MslWriter.Write(path, new[] { Entry("PEPTIDEK", 450.0, neutralLoss: -18.010565) });

		Assert.That(MslReader.ReadHeaderOnly(path).FileFlags & MslFormat.FileFlagHasExtAnnotations, Is.Not.EqualTo(0));
		using MslLibrary library = MslLibrary.Load(path);
		MslFragmentIon read = library.GetEntry(0, includeFragments: true)!.MatchedFragmentIons.Single();
		Assert.That(read.NeutralLoss, Is.EqualTo(-18.010565));
		Assert.That(read.Annotation, Is.EqualTo("y4+18.0106"));
	}

	// ── Merge: the winner follows the policy, never the m/z heap ─────────────

	/// <summary>
	/// KeepFirst keeps the first file's copy when the copies tie on m/z (every reversed decoy
	/// ties its target). Before the fix the heap's order chose, and B's copy won.
	/// </summary>
	[Test]
	public void Merge_KeepFirst_TiedMz_KeepsTheFirstFilesCopy()
	{
		string a = PathFor("tie_a"), b = PathFor("tie_b"), output = PathFor("tie_out");
		MslWriter.Write(a, new[] { Entry("PEPTIDEP", 300.0, rt: 1f), Entry("PEPTIDEK", 300.0, rt: 1f) });
		MslWriter.Write(b, new[] { Entry("PEPTIDEK", 300.0, rt: 99f) });

		MslMergeResult result = MslMerger.Merge(new[] { a, b }, output, MslMergeConflictPolicy.KeepFirst);

		using MslLibrary merged = MslLibrary.Load(output);
		Assert.That(merged.TryGetEntry("PEPTIDEK", 2, out MslLibraryEntry? kept), Is.True);
		Assert.That(kept!.RetentionTime, Is.EqualTo(1f));
		Assert.That(result.OutputEntryCount, Is.EqualTo(2));
		Assert.That(result.DuplicatesSkipped, Is.EqualTo(1));
	}

	/// <summary>
	/// KeepFirst keeps the first file's copy even when the later file's copy has the lower m/z.
	/// Before the fix the lower m/z was emitted first and won.
	/// </summary>
	[Test]
	public void Merge_KeepFirst_CopiesApartInMz_KeepsTheFirstFilesCopy()
	{
		string a = PathFor("first_a"), b = PathFor("first_b"), output = PathFor("first_out");
		MslWriter.Write(a, new[] { Entry("PEPTIDEK", 300.002, rt: 1f) });
		MslWriter.Write(b, new[] { Entry("PEPTIDEK", 300.000, rt: 99f) });

		MslMerger.Merge(new[] { a, b }, output, MslMergeConflictPolicy.KeepFirst);

		using MslLibrary merged = MslLibrary.Load(output);
		Assert.That(merged.PrecursorCount, Is.EqualTo(1));
		Assert.That(merged.GetEntry(0, includeFragments: false)!.RetentionTime, Is.EqualTo(1f));
	}

	/// <summary>
	/// KeepLowestQValue resolves two copies of one Name whose m/z differ by 0.002 Da.
	/// Before the fix it compared only copies within 0.001 Da, and kept both.
	/// </summary>
	[Test]
	public void Merge_KeepLowestQValue_CopiesApartInMz_KeepsOnlyTheBestCopy()
	{
		string a = PathFor("q_a"), b = PathFor("q_b"), output = PathFor("q_out");
		MslWriter.Write(a, new[] { Entry("PEPTIDEK", 300.000, qValue: 0.05f) });
		MslWriter.Write(b, new[] { Entry("PEPTIDEK", 300.002, qValue: 0.01f) });

		MslMergeResult result = MslMerger.Merge(new[] { a, b }, output, MslMergeConflictPolicy.KeepLowestQValue);

		Assert.That(result.OutputEntryCount, Is.EqualTo(1));
		Assert.That(result.DuplicatesSkipped, Is.EqualTo(1));
		using MslLibrary merged = MslLibrary.Load(output);
		Assert.That(merged.GetEntry(0, includeFragments: false)!.QValue, Is.EqualTo(0.01f));
	}

	/// <summary>KeepLast keeps the last file's copy whatever the m/z order (guards the rewrite).</summary>
	[Test]
	public void Merge_KeepLast_CopiesApartInMz_KeepsTheLastFilesCopy()
	{
		string a = PathFor("last_a"), b = PathFor("last_b"), output = PathFor("last_out");
		MslWriter.Write(a, new[] { Entry("PEPTIDEK", 300.000, rt: 1f) });
		MslWriter.Write(b, new[] { Entry("PEPTIDEK", 299.998, rt: 99f) });

		MslMerger.Merge(new[] { a, b }, output, MslMergeConflictPolicy.KeepLast);

		using MslLibrary merged = MslLibrary.Load(output);
		Assert.That(merged.PrecursorCount, Is.EqualTo(1));
		Assert.That(merged.GetEntry(0, includeFragments: false)!.RetentionTime, Is.EqualTo(99f));
	}

	/// <summary>
	/// KeepLast writes a Name once even when the last file holds it twice. Before the fix every
	/// copy in the last file was written, so the deduplicated output held the Name twice.
	/// </summary>
	[Test]
	public void Merge_KeepLast_ANameTheLastFileHoldsTwice_IsWrittenOnce()
	{
		string a = PathFor("twice_a"), b = PathFor("twice_b"), output = PathFor("twice_out");
		MslWriter.Write(a, new[] { Entry("PEPTIDEK", 300.0, rt: 1f) });
		MslWriter.Write(b, new[] { Entry("PEPTIDEK", 300.0, rt: 50f), Entry("PEPTIDEK", 300.001, rt: 99f) });

		MslMergeResult result = MslMerger.Merge(new[] { a, b }, output, MslMergeConflictPolicy.KeepLast);

		Assert.That(result.OutputEntryCount, Is.EqualTo(1));
		Assert.That(result.DuplicatesSkipped, Is.EqualTo(2));
		using MslLibrary merged = MslLibrary.Load(output);
		Assert.That(merged.PrecursorCount, Is.EqualTo(1));
		Assert.That(merged.GetEntry(0, includeFragments: false)!.RetentionTime, Is.EqualTo(99f),
			"The last copy of the last file wins.");
	}

	/// <summary>
	/// An output path that is one of the sources is refused before any work, and the source is
	/// untouched. Before the fix Windows failed late (the source is held open), and elsewhere
	/// the output replaced the source.
	/// </summary>
	[Test]
	public void Merge_AnOutputThatIsASource_IsRefusedAndTheSourceUntouched()
	{
		string a = PathFor("self_a"), b = PathFor("self_b");
		MslWriter.Write(a, new[] { Entry("PEPTIDEK", 300.0) });
		MslWriter.Write(b, new[] { Entry("PEPTIDER", 310.0) });
		byte[] before = File.ReadAllBytes(a);

		Assert.That(() => MslMerger.Merge(new[] { a, b }, a), Throws.ArgumentException);
		Assert.That(File.ReadAllBytes(a), Is.EqualTo(before));
	}

	// ── Compression level: one rule for every writer ─────────────────────────

	/// <summary>
	/// WriteStreaming and Merge refuse a level outside [0, 22], as Write always did.
	/// Before the fix -1 wrote uncompressed and 23 wrote compressed.
	/// </summary>
	[TestCase(-1)]
	[TestCase(23)]
	public void WriteStreaming_AndMerge_RefuseALevelOutside0To22(int level)
	{
		string source = PathFor($"level_src_{level}");
		string path = PathFor($"level_{level}");
		MslWriter.Write(source, new[] { Entry("PEPTIDEK", 300.0) });

		Assert.That(() => MslWriter.WriteStreaming(path, new[] { Entry("PEPTIDEK", 300.0) }, level),
			Throws.TypeOf<ArgumentOutOfRangeException>());
		Assert.That(() => MslMerger.Merge(new[] { source }, path, compressionLevel: level),
			Throws.TypeOf<ArgumentOutOfRangeException>());
		Assert.That(File.Exists(path), Is.False);
	}

	// ── The window result's Dispose ──────────────────────────────────────────

	/// <summary>
	/// Disposing a window result more than once, or disposing a copy too, returns its pooled
	/// array once. Before the fix each call returned it to ArrayPool.Shared, and the next two
	/// queries were handed the same array, so one overwrote the other's results.
	/// </summary>
	[Test]
	public void QueryWindow_DisposedTwice_TheNextTwoResultsDoNotShareAnArray()
	{
		string path = PathFor(nameof(QueryWindow_DisposedTwice_TheNextTwoResultsDoNotShareAnArray));
		MslWriter.Write(path, Enumerable.Range(0, 40).Select(i => Entry($"PEPTIDE{i:D2}K", 400.0 + i)).ToList());
		using MslLibrary library = MslLibrary.Load(path);

		MslWindowResults first = library.QueryWindow(399f, 500f, 0f, 100f);
		Assert.That(first.Count, Is.EqualTo(40));
		MslWindowResults copy = first;
		first.Dispose();
		copy.Dispose();
		first.Dispose();
		Assert.That(first.Entries.Length, Is.EqualTo(0), "A disposed result shows no entries.");

		using MslWindowResults second = library.QueryWindow(399f, 500f, 0f, 100f);
		using MslWindowResults third = library.QueryWindow(399f, 500f, 0f, 100f);
		Assert.That(Unsafe.AreSame(
				ref MemoryMarshal.GetReference(second.Entries),
				ref MemoryMarshal.GetReference(third.Entries)),
			Is.False, "Two live results must not share one pooled array.");
	}

	// ── Unconvertible spectra ────────────────────────────────────────────────

	/// <summary>
	/// A spectrum FromLibrarySpectrum cannot convert (here charge 0) is an ArgumentException that
	/// names it, and nothing is written. Before the fix it was a NullReferenceException.
	/// </summary>
	[Test]
	public void AnUnconvertibleSpectrum_IsAnArgumentExceptionNamingIt()
	{
		string path = PathFor(nameof(AnUnconvertibleSpectrum_IsAnArgumentExceptionNamingIt));
		LibrarySpectrum good = Spectrum("PEPTIDEK", 2);
		LibrarySpectrum bad = Spectrum("PEPTIDER", 0);

		Assert.That(() => MslWriter.WriteFromLibrarySpectra(path, new[] { good, bad }),
			Throws.ArgumentException.With.Message.Contains("spectra[1]"));
		Assert.That(() => MslLibrary.SaveFromLibrarySpectra(path, new[] { good, bad }), Throws.ArgumentException);
		Assert.That(() => MslConverter.FromMspLibrarySpectra(new[] { good, bad }), Throws.ArgumentException);
		Assert.That(File.Exists(path), Is.False);

		string source = PathFor("update_source");
		MslWriter.WriteFromLibrarySpectra(source, new[] { good });
		using MslLibrary library = MslLibrary.Load(source);
		Assert.That(() => library.UpdateAndSave(new[] { bad }, path),
			Throws.ArgumentException.With.Message.Contains("newSpectra[0]"));
		Assert.That(File.Exists(path), Is.False);
	}

	// ── Values the writer used to store wrongly ──────────────────────────────

	/// <summary>
	/// A precursor m/z that is not a finite float is refused. Before the fix it was written,
	/// and then no query or GetAllEntries could ever return the entry.
	/// </summary>
	[TestCase(double.NaN)]
	[TestCase(double.PositiveInfinity)]
	[TestCase(1e40)] // finite as a double, infinite as the stored float
	public void APrecursorMzThatIsNotAFiniteFloat_IsRefused(double mz)
	{
		string path = PathFor($"mz_{mz}");
		Assert.That(() => MslWriter.Write(path, new[] { Entry("PEPTIDEK", mz) }), Throws.ArgumentException);
		Assert.That(() => MslWriter.WriteStreaming(path, new[] { Entry("PEPTIDEK", mz) }), Throws.ArgumentException);
		Assert.That(File.Exists(path), Is.False);
	}

	/// <summary>
	/// A fragment charge outside the stored byte is refused. Before the fix -1 was stored as 255.
	/// </summary>
	[TestCase(-1)]
	[TestCase(256)]
	public void AFragmentChargeOutside0To255_IsRefused(int charge)
	{
		string path = PathFor($"fragcharge_{charge}");
		Assert.That(() => MslWriter.Write(path, new[] { Entry("PEPTIDEK", 300.0, fragmentCharge: charge) }),
			Throws.ArgumentException);
		Assert.That(() => MslWriter.WriteStreaming(path, new[] { Entry("PEPTIDEK", 300.0, fragmentCharge: charge) }),
			Throws.ArgumentException);
		Assert.That(File.Exists(path), Is.False);
	}

	// ── Temp-file names ──────────────────────────────────────────────────────

	/// <summary>
	/// A write's temp files carry a per-call token, so a write does not touch another writer's
	/// temp file. Before the fix every write used "&lt;out&gt;.tmp~": on Windows it failed while
	/// another writer held that file, and elsewhere it truncated and then renamed it away.
	/// </summary>
	[Test]
	public void Write_DoesNotTouchAnotherWritersTempFile()
	{
		string path = PathFor(nameof(Write_DoesNotTouchAnotherWritersTempFile));
		string otherWritersTemp = path + ".tmp~";
		using (var held = new FileStream(otherWritersTemp, FileMode.Create, FileAccess.Write, FileShare.None))
		{
			held.WriteByte(42);
			Assert.That(() => MslWriter.Write(path, new[] { Entry("PEPTIDEK", 300.0) }), Throws.Nothing);
			Assert.That(() => MslWriter.WriteStreaming(path, new[] { Entry("PEPTIDEK", 300.0) }), Throws.Nothing);
		}

		Assert.That(File.ReadAllBytes(otherWritersTemp), Is.EqualTo(new byte[] { 42 }));
		Assert.That(Directory.GetFiles(_dir, Path.GetFileName(path) + ".*~"), Is.EquivalentTo(new[] { otherWritersTemp }),
			"The writes left no temp files of their own.");
		using MslLibrary written = MslLibrary.Load(path);
		Assert.That(written.PrecursorCount, Is.EqualTo(1));
	}

	// ── Disposal during enumeration ──────────────────────────────────────────

	/// <summary>
	/// Disposing the library while GetAllEntries is being enumerated ends the enumeration, as
	/// its doc says. Before the fix the next MoveNext threw NullReferenceException.
	/// </summary>
	[TestCase(false)]
	[TestCase(true)]
	public void GetAllEntries_DisposedMidEnumeration_StopsCleanly(bool indexOnly)
	{
		string path = PathFor($"enumerate_{indexOnly}");
		MslWriter.Write(path, new[] { Entry("PEPTIDEK", 300.0), Entry("PEPTIDER", 310.0), Entry("PEPTIDES", 320.0) });
		MslLibrary library = indexOnly ? MslLibrary.LoadIndexOnly(path) : MslLibrary.Load(path);

		using IEnumerator<MslLibraryEntry> entries = library.GetAllEntries().GetEnumerator();
		Assert.That(entries.MoveNext(), Is.True);
		library.Dispose();

		bool more = true;
		Assert.That(() => more = entries.MoveNext(), Throws.Nothing);
		Assert.That(more, Is.False);
	}

	// ── Calibration slope ────────────────────────────────────────────────────

	/// <summary>
	/// A slope that cannot calibrate is refused. Before the fix a NaN slope made every index RT
	/// NaN, so every RT window matched, and a zero slope put every RT on the intercept.
	/// </summary>
	[TestCase(double.NaN)]
	[TestCase(double.PositiveInfinity)]
	[TestCase(0.0)]
	public void WithCalibratedRetentionTimes_ASlopeThatCannotCalibrate_IsRefused(double slope)
	{
		string path = PathFor($"slope_{slope}");
		MslWriter.Write(path, new[] { Entry("PEPTIDEK", 300.0) });
		using MslLibrary library = MslLibrary.Load(path);

		Assert.That(() => library.WithCalibratedRetentionTimes(slope, 1.0), Throws.TypeOf<ArgumentOutOfRangeException>());
	}

	// ── Header version ───────────────────────────────────────────────────────

	/// <summary>
	/// ReadHeaderOnly refuses a format version this reader does not know, as a full load does.
	/// Before the fix it returned the header of a version-99 file as if it were fine.
	/// </summary>
	[Test]
	public void ReadHeaderOnly_AnUnknownVersion_IsAFormatException()
	{
		string path = PathFor(nameof(ReadHeaderOnly_AnUnknownVersion_IsAFormatException));
		MslWriter.Write(path, new[] { Entry("PEPTIDEK", 300.0) });
		byte[] bytes = File.ReadAllBytes(path);
		BitConverter.GetBytes(99).CopyTo(bytes, 4); // MslFileHeader.FormatVersion, offset 4
		File.WriteAllBytes(path, bytes);

		Assert.That(() => MslReader.ReadHeaderOnly(path), Throws.TypeOf<FormatException>());
	}

	// ── Proton mass ──────────────────────────────────────────────────────────

	/// <summary>The proteoform code uses mzLib's one proton constant, not a rounded copy.</summary>
	[Test]
	public void ProteoformNeutralMass_UsesMzLibsProtonConstant()
	{
		Assert.That(MslProteoformIndex.ProtonMass, Is.EqualTo(Constants.ProtonMass));
	}

	// ── A missing retention time is NaN, never iRT 0 ─────────────────────────
	// 0.0 is a real iRT value (where GAGSSEPVTGLDAK elutes), so it cannot also mean "none".

	/// <summary>An entry built without a time has none. Before the fix it defaulted to iRT 0.</summary>
	[Test]
	public void ADefaultEntry_HasNoRetentionTime()
	{
		Assert.That(new MslLibraryEntry().RetentionTime, Is.NaN);
	}

	/// <summary>
	/// A spectrum without a time (an MSP entry with no RT=) gives NaN through every conversion and
	/// round-trips as NaN; back out as a LibrarySpectrum it is null, so an MSP export omits RT=.
	/// Before the fix every such spectrum became iRT 0, and RT=0 was then written out as a real time.
	/// </summary>
	[Test]
	public void ASpectrumWithoutATime_IsNaNThroughEveryConversion_AndNullBackOut()
	{
		LibrarySpectrum noTime = Spectrum("PEPTIDEK", 2, rt: null);

		string viaWriter = PathFor("nort_writer");
		MslWriter.WriteFromLibrarySpectra(viaWriter, new[] { noTime });
		string viaSave = PathFor("nort_save");
		MslLibrary.SaveFromLibrarySpectra(viaSave, new[] { noTime });
		string msp = Path.Combine(_dir, "nort.msp");
		File.WriteAllText(msp, noTime.ToString());
		Assert.That(File.ReadAllText(msp), Does.Not.Contain("RT="));

		foreach (string path in new[] { viaWriter, viaSave })
		{
			using MslLibrary library = MslLibrary.Load(path);
			MslLibraryEntry entry = library.GetEntry(0, includeFragments: true)!;
			Assert.That(entry.RetentionTime, Is.NaN, path);
			Assert.That(entry.ToLibrarySpectrum().RetentionTime, Is.Null, path);
		}

		using MslLibrary converted = MslConverter.FromMspFile(msp);
		MslLibraryEntry fromMsp = converted.GetEntry(0, includeFragments: true)!;
		Assert.That(fromMsp.RetentionTime, Is.NaN);
		Assert.That(fromMsp.ToLibrarySpectrum().ToString(), Does.Not.Contain("RT="),
			"Exporting the entry again must not invent RT=0 or RT=NaN.");
	}

	/// <summary>
	/// An RT window keeps an entry with no time whatever its bounds, and keeps an entry at iRT 0
	/// only when the window contains 0. Before the fix the "no time" entry was an iRT-0 entry,
	/// hidden from every window but one around 0.
	/// </summary>
	[Test]
	public void QueryWindow_AnEntryWithNoTimePassesAnyWindow_AnEntryAtZeroOnlyWindowsWithZero()
	{
		string path = PathFor(nameof(QueryWindow_AnEntryWithNoTimePassesAnyWindow_AnEntryAtZeroOnlyWindowsWithZero));
		MslWriter.Write(path, new[] { Entry("PEPTIDEN", 300.0, rt: float.NaN), Entry("PEPTIDEZ", 301.0, rt: 0f) });
		using MslLibrary library = MslLibrary.Load(path);

		string[] Names(float rtLow, float rtHigh)
		{
			using MslWindowResults window = library.QueryWindow(299f, 302f, rtLow, rtHigh);
			var names = new List<string>();
			foreach (MslPrecursorIndexEntry e in window.Entries)
				names.Add(library.GetEntry(e.PrecursorIdx, includeFragments: false)!.FullSequence);
			return names.ToArray();
		}

		Assert.That(Names(50f, 60f), Is.EquivalentTo(new[] { "PEPTIDEN" }));
		Assert.That(Names(-5f, 5f), Is.EquivalentTo(new[] { "PEPTIDEN", "PEPTIDEZ" }));
	}

	/// <summary>Calibration keeps a missing time missing, and the iRT range ignores it.</summary>
	[Test]
	public void ANaNTime_StaysNaNUnderCalibration_AndIsLeftOutOfTheIrtRange()
	{
		string path = PathFor(nameof(ANaNTime_StaysNaNUnderCalibration_AndIsLeftOutOfTheIrtRange));
		MslWriter.Write(path, new[]
		{
			Entry("PEPTIDEN", 300.0, rt: float.NaN), Entry("PEPTIDEA", 301.0, rt: 10f), Entry("PEPTIDEB", 302.0, rt: 20f)
		});
		using MslLibrary library = MslLibrary.Load(path);

		MslIndexStatistics stats = library.IndexStatisticsForTesting;
		Assert.That(stats.MinIrt, Is.EqualTo(10.0));
		Assert.That(stats.MaxIrt, Is.EqualTo(20.0));

		using MslLibrary calibrated = library.WithCalibratedRetentionTimes(2.0, 5.0);
		float[] irts = calibrated.QueryMzWindow(299f, 303f).ToArray().Select(e => e.Irt).ToArray();
		Assert.That(irts[0], Is.NaN);
		Assert.That(irts.Skip(1), Is.EqualTo(new[] { 25f, 45f }));
	}

	/// <summary>
	/// The writer refuses an infinite retention time (NaN alone means "no time"), and
	/// ValidateEntries reports it.
	/// </summary>
	[TestCase(double.PositiveInfinity)]
	[TestCase(double.NegativeInfinity)]
	[TestCase(1e40)] // finite as a double, infinite as the stored float
	public void AnInfiniteRetentionTime_IsRefused(double rt)
	{
		string path = PathFor($"infrt_{rt}");
		var entry = Entry("PEPTIDEK", 300.0);
		entry.RetentionTime = rt;

		Assert.That(() => MslWriter.Write(path, new[] { entry }), Throws.ArgumentException);
		Assert.That(() => MslWriter.WriteStreaming(path, new[] { entry }), Throws.ArgumentException);
		Assert.That(File.Exists(path), Is.False);
		Assert.That(MslWriter.ValidateEntries(new[] { entry }), Has.Some.Contains("RetentionTime"));
	}

	/// <summary>A NaN retention time is written and read back as NaN.</summary>
	[Test]
	public void ANaNRetentionTime_IsWrittenAndReadBack()
	{
		string path = PathFor(nameof(ANaNRetentionTime_IsWrittenAndReadBack));
		MslWriter.WriteStreaming(path, new[] { Entry("PEPTIDEK", 300.0, rt: float.NaN) });
		using MslLibrary library = MslLibrary.LoadIndexOnly(path);
		Assert.That(library.GetEntry(0, includeFragments: false)!.RetentionTime, Is.NaN);
	}

	// ── UpdateAndSave and missing times ──────────────────────────────────────

	/// <summary>A library of five entries at iRT 10, 20, ... 50, one fragment each.</summary>
	private static MslLibrary AnchorLibrary(string name, params MslLibraryEntry[] extra)
	{
		string path = PathFor(name);
		var entries = Enumerable.Range(1, 5).Select(i => Entry($"ANCHOR{i}K", 400.0 + i, rt: 10f * i)).ToList();
		entries.AddRange(extra);
		MslWriter.Write(path, entries);
		return MslLibrary.Load(path);
	}

	/// <summary>Observed spectra for the five anchors, at 2, 4, ... 10 minutes (iRT = 5 × minutes).</summary>
	private static List<LibrarySpectrum> AnchorSpectra() =>
		Enumerable.Range(1, 5).Select(i => Spectrum($"ANCHOR{i}K", 2, rt: 2.0 * i, precursorMz: 400.0 + i)).ToList();

	/// <summary>
	/// Only pairs of real times are anchors: an observed time that is missing, 0 or negative (psmtsv
	/// writes -1), or a library time that is NaN, does not pull the line. Before the fix a missing
	/// observed time became an anchor at 0 minutes, and a NaN library time made the fit NaN.
	/// </summary>
	[Test]
	public void UpdateAndSave_OnlyRealTimesAreAnchors()
	{
		using MslLibrary clean = AnchorLibrary("anchors_clean");
		MslUpdateResult expected = clean.UpdateAndSave(AnchorSpectra(), PathFor("anchors_clean_out"));
		Assert.That(expected.RtNormalisationApplied, Is.True);
		Assert.That(expected.RtSlope, Is.EqualTo(5.0).Within(1e-4));
		Assert.That(expected.RtIntercept, Is.EqualTo(0.0).Within(1e-3));

		using MslLibrary noisy = AnchorLibrary("anchors_noisy",
			Entry("MISSINGK", 410.0, rt: 30f), Entry("ZEROK", 411.0, rt: 30f), Entry("NEGATIVEK", 412.0, rt: 30f),
			Entry("NANLIBK", 413.0, rt: float.NaN));
		List<LibrarySpectrum> spectra = AnchorSpectra();
		spectra.Add(Spectrum("MISSINGK", 2, rt: null, precursorMz: 410.0));
		spectra.Add(Spectrum("ZEROK", 2, rt: 0.0, precursorMz: 411.0));
		spectra.Add(Spectrum("NEGATIVEK", 2, rt: -1.0, precursorMz: 412.0));
		spectra.Add(Spectrum("NANLIBK", 2, rt: 7.0, precursorMz: 413.0));

		MslUpdateResult result = noisy.UpdateAndSave(spectra, PathFor("anchors_noisy_out"));

		Assert.That(result.AnchorCount, Is.EqualTo(5));
		Assert.That(result.RtNormalisationApplied, Is.True);
		Assert.That(result.RtSlope, Is.EqualTo(expected.RtSlope).Within(1e-9));
		Assert.That(result.RtIntercept, Is.EqualTo(expected.RtIntercept).Within(1e-9));
	}

	/// <summary>
	/// A better spectrum without a time keeps the replaced entry's time and flag, and keeps its
	/// protein and gene names. Before the fix the time became NormaliseRt(0), an invented value,
	/// and the names were dropped.
	/// </summary>
	[Test]
	public void UpdateAndSave_AReplacementWithoutATime_KeepsTheKnownTime_AndTheNames()
	{
		var original = Entry("KEPTK", 410.0, rt: 42f);
		original.ProteinAccession = "P12345";
		original.ProteinName = "Kept protein";
		original.GeneName = "KEPT1";
		using MslLibrary library = AnchorLibrary("replace_keep", original);

		List<LibrarySpectrum> spectra = AnchorSpectra();
		spectra.Add(Spectrum("KEPTK", 2, rt: null, fragments: 3, precursorMz: 410.0));
		string output = PathFor("replace_keep_out");
		MslUpdateResult result = library.UpdateAndSave(spectra, output);
		Assert.That(result.ReplacedCount, Is.EqualTo(1));

		using MslLibrary updated = MslLibrary.Load(output);
		Assert.That(updated.TryGetEntry("KEPTK", 2, out MslLibraryEntry? kept), Is.True);
		Assert.That(kept!.MatchedFragmentIons, Has.Count.EqualTo(3), "The better spectrum replaced the entry.");
		Assert.That(kept.RetentionTime, Is.EqualTo(42.0));
		Assert.That(kept.RtIsCalibrated, Is.False);
		Assert.That(kept.ProteinAccession, Is.EqualTo("P12345"));
		Assert.That(kept.ProteinName, Is.EqualTo("Kept protein"));
		Assert.That(kept.GeneName, Is.EqualTo("KEPT1"));
		Assert.That(kept.QValue, Is.NaN, "The replaced spectrum's q-value does not describe the new one.");
	}

	/// <summary>A novel entry without a time has none. Before the fix it was stored at NormaliseRt(0).</summary>
	[Test]
	public void UpdateAndSave_ANovelEntryWithoutATime_HasNone()
	{
		using MslLibrary library = AnchorLibrary("novel_nan");
		List<LibrarySpectrum> spectra = AnchorSpectra();
		spectra.Add(Spectrum("NOVELK", 2, rt: null, precursorMz: 420.0));
		string output = PathFor("novel_nan_out");
		MslUpdateResult result = library.UpdateAndSave(spectra, output);
		Assert.That(result.AddedCount, Is.EqualTo(1));

		using MslLibrary updated = MslLibrary.Load(output);
		Assert.That(updated.TryGetEntry("NOVELK", 2, out MslLibraryEntry? novel), Is.True);
		Assert.That(novel!.RetentionTime, Is.NaN);
	}

	/// <summary>
	/// A fit that is not finite (every anchor at one observed time) is not reported as applied.
	/// Before the fix RtNormalisationApplied was true with a NaN slope.
	/// </summary>
	[Test]
	public void UpdateAndSave_ADegenerateFit_IsNotApplied()
	{
		using MslLibrary library = AnchorLibrary("degenerate");
		List<LibrarySpectrum> spectra =
			Enumerable.Range(1, 5).Select(i => Spectrum($"ANCHOR{i}K", 2, rt: 6.0, precursorMz: 400.0 + i)).ToList();

		MslUpdateResult result = library.UpdateAndSave(spectra, PathFor("degenerate_out"));

		Assert.That(result.RtNormalisationApplied, Is.False);
		Assert.That(double.IsFinite(result.RtSlope) && double.IsFinite(result.RtIntercept), Is.True);
	}

	/// <summary>
	/// PredictFragments refuses an index-only library, where eviction would discard its work.
	/// Before the fix it reported success and kept nothing.
	/// </summary>
	[Test]
	public void PredictFragments_RefusesAnIndexOnlyLibrary()
	{
		string path = PathFor(nameof(PredictFragments_RefusesAnIndexOnlyLibrary));
		var predicted = Entry("PREDICTK", 300.0);
		predicted.Source = MslFormat.SourceType.Predicted;
		predicted.MatchedFragmentIons.Clear();
		MslWriter.Write(path, new[] { predicted });

		using MslLibrary indexOnly = MslLibrary.LoadIndexOnly(path);
		Assert.That(() => indexOnly.PredictFragments(_ => new Dictionary<string, LibrarySpectrum>()),
			Throws.InvalidOperationException);
	}

	// ── IrtCalibrationModel and missing times ────────────────────────────────

	/// <summary>
	/// A missing time maps to a missing time, and is never inside the anchors' span. Before the
	/// fix ToIrt(NaN) and ToRtMinutes(NaN) threw IndexOutOfRangeException and IsExtrapolated(NaN)
	/// was false.
	/// </summary>
	[Test]
	public void IrtCalibrationModel_ANaNTime_MapsToNaN_AndIsExtrapolated()
	{
		IrtCalibrationModel model = IrtCalibration.Line((new RtMinutes(10), new Irt(0)), (new RtMinutes(20), new Irt(100)));

		Assert.That(model.ToIrt(new RtMinutes(double.NaN)).Value, Is.NaN);
		Assert.That(model.ToRtMinutes(new Irt(double.NaN)).Value, Is.NaN);
		Assert.That(model.IsExtrapolated(new RtMinutes(double.NaN)), Is.True);
		Assert.That(model.ToIrt(new RtMinutes(15)).Value, Is.EqualTo(50.0).Within(1e-9), "A real time still maps.");
	}
}
