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
/// Tests for the retention-time status (PrecursorFlags bits 4-5) and the struck-by-prediction
/// flag (bit 6). Like the entrapment bit, they were added without a format-version bump: older
/// readers ignore them. Before this change every reader dropped them and every writer stored 0,
/// so any read and rewrite (a merge, UpdateAndSave, Save of GetAllEntries) erased them.
/// </summary>
[TestFixture]
public sealed class TestMslRetentionTimeStatus
{
	private static readonly MslFormat.RetentionTimeStatus[] AllStatuses =
		Enum.GetValues<MslFormat.RetentionTimeStatus>();

	private string _dir = null!;

	[OneTimeSetUp]
	public void OneTimeSetUp()
	{
		_dir = Path.Combine(Path.GetTempPath(), "MslRetentionTimeStatus_" + Guid.NewGuid().ToString("N"));
		Directory.CreateDirectory(_dir);
	}

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		try { Directory.Delete(_dir, recursive: true); } catch (IOException) { }
	}

	private string TempPath(string name) => Path.Combine(_dir, name + ".msl");

	private static MslLibraryEntry Entry(string seq, int charge, double mz, double rt,
		MslFormat.RetentionTimeStatus status = MslFormat.RetentionTimeStatus.None, bool struck = false,
		bool isDecoy = false, bool isEntrapment = false, bool rtIsCalibrated = false) => new()
	{
		FullSequence = seq,
		BaseSequence = seq,
		PrecursorMz = mz,
		ChargeState = charge,
		RetentionTime = rt,
		RtStatus = status,
		RtStruckByPrediction = struck,
		IsDecoy = isDecoy,
		IsEntrapment = isEntrapment,
		RtIsCalibrated = rtIsCalibrated,
		ProteinAccession = "P" + seq,
		DissociationType = DissociationType.HCD,
		Nce = 28,
		MatchedFragmentIons = new List<MslFragmentIon>
		{
			new() { Mz = 200f + seq.Length, Intensity = 1.0f, ProductType = ProductType.b, FragmentNumber = 2, Charge = 1 },
			new() { Mz = 300f + seq.Length, Intensity = 0.5f, ProductType = ProductType.y, FragmentNumber = 3, Charge = 1 },
		}
	};

	/// <summary>
	/// One entry per status x struck combination (8), with the other flags varied alongside, so
	/// a bit that leaks into its neighbour shows up.
	/// </summary>
	private static List<MslLibraryEntry> EveryCombination()
	{
		var entries = new List<MslLibraryEntry>();
		int i = 0;
		foreach (MslFormat.RetentionTimeStatus status in AllStatuses)
			foreach (bool struck in new[] { false, true })
			{
				string seq = "PEPTIDE" + (char)('A' + i) + "K";
				entries.Add(Entry(seq, 2, 400.0 + 11.0 * i, 10.0 + i, status, struck,
					isDecoy: i % 2 == 1, isEntrapment: i % 3 == 0, rtIsCalibrated: i % 4 == 2));
				i++;
			}
		return entries;
	}

	private static void WriteWith(string writer, string path, List<MslLibraryEntry> entries, int compressionLevel)
	{
		if (writer == "Write")
			MslWriter.Write(path, entries, compressionLevel);
		else
			MslWriter.WriteStreaming(path, entries, compressionLevel);
	}

	private static void AssertSameFlags(MslLibraryEntry actual, MslLibraryEntry expected, string where)
	{
		Assert.That(actual.RtStatus, Is.EqualTo(expected.RtStatus), $"{where}: RtStatus");
		Assert.That(actual.RtStruckByPrediction, Is.EqualTo(expected.RtStruckByPrediction), $"{where}: RtStruckByPrediction");
		Assert.That(actual.IsDecoy, Is.EqualTo(expected.IsDecoy), $"{where}: IsDecoy");
		Assert.That(actual.IsEntrapment, Is.EqualTo(expected.IsEntrapment), $"{where}: IsEntrapment");
		Assert.That(actual.RtIsCalibrated, Is.EqualTo(expected.RtIsCalibrated), $"{where}: RtIsCalibrated");
	}

	// ── The flags byte ────────────────────────────────────────────────────────

	[Test]
	public void EveryStatusAndStrike_RoundTripsThroughTheFlagsByte_WithEveryOtherFlag()
	{
		bool[] tf = { false, true };
		int combinations = 0;
		foreach (MslFormat.RetentionTimeStatus status in AllStatuses)
		foreach (bool struck in tf)
		foreach (bool decoy in tf)
		foreach (bool proteotypic in tf)
		foreach (bool calibrated in tf)
		foreach (bool entrapment in tf)
		{
			byte flags = MslFormat.EncodePrecursorFlags(decoy, proteotypic, calibrated, entrapment, status, struck);

			Assert.That(MslFormat.DecodeRetentionTimeStatus(flags), Is.EqualTo(status));
			Assert.That(MslFormat.DecodeRtStruckByPrediction(flags), Is.EqualTo(struck));
			Assert.That(MslFormat.DecodePrecursorFlags(flags), Is.EqualTo((decoy, proteotypic, calibrated)));
			Assert.That(MslFormat.DecodeIsEntrapment(flags), Is.EqualTo(entrapment));
			Assert.That(flags & 0b_1000_0000, Is.Zero, "bit 7 stays reserved");
			combinations++;
		}
		Assert.That(combinations, Is.EqualTo(4 * 2 * 16));
	}

	[Test]
	public void OlderEncodeOverloads_WriteNoStatus()
	{
		byte three = MslFormat.EncodePrecursorFlags(true, true, true);
		byte four = MslFormat.EncodePrecursorFlags(true, true, true, true);

		Assert.That(three & 0b_1111_0000, Is.Zero);
		Assert.That(four & 0b_1111_0000, Is.Zero);
		Assert.That(MslFormat.DecodeRetentionTimeStatus(four), Is.EqualTo(MslFormat.RetentionTimeStatus.None));
		Assert.That(MslFormat.DecodeRtStruckByPrediction(four), Is.False);
	}

	[Test]
	public void ReservedBit7_IsIgnoredOnRead()
	{
		byte flags = (byte)(0b_1000_0000 | MslFormat.EncodePrecursorFlags(false, false, false, false,
			MslFormat.RetentionTimeStatus.Confirmed, rtStruckByPrediction: false));

		Assert.That(MslFormat.DecodeRetentionTimeStatus(flags), Is.EqualTo(MslFormat.RetentionTimeStatus.Confirmed));
		Assert.That(MslFormat.DecodeRtStruckByPrediction(flags), Is.False);
	}

	[Test]
	public void Encoder_RefusesAnUndefinedStatus() =>
		Assert.That(() => MslFormat.EncodePrecursorFlags(false, false, false, false,
				(MslFormat.RetentionTimeStatus)4, rtStruckByPrediction: false),
			Throws.TypeOf<ArgumentOutOfRangeException>().With.Property("ParamName").EqualTo("rtStatus"));

	// ── Files ─────────────────────────────────────────────────────────────────

	[TestCase("Write", 0, false)]
	[TestCase("Write", 0, true)]
	[TestCase("Write", 3, false)]
	[TestCase("Write", 3, true)]
	[TestCase("WriteStreaming", 0, false)]
	[TestCase("WriteStreaming", 0, true)]
	[TestCase("WriteStreaming", 3, false)]
	[TestCase("WriteStreaming", 3, true)]
	public void StatusAndStrike_SurviveWriteAndRead(string writer, int compressionLevel, bool indexOnly)
	{
		List<MslLibraryEntry> expected = EveryCombination();
		string path = TempPath($"roundtrip_{writer}_{compressionLevel}_{indexOnly}");
		WriteWith(writer, path, EveryCombination(), compressionLevel);

		using MslLibrary lib = indexOnly ? MslLibrary.LoadIndexOnly(path) : MslLibrary.Load(path);
		foreach (MslLibraryEntry e in expected)
		{
			Assert.That(lib.TryGetEntry(e.FullSequence, e.ChargeState, out MslLibraryEntry? read), Is.True, e.FullSequence);
			AssertSameFlags(read!, e, e.FullSequence);
		}

		// The metadata-only read (index-only skeletons) carries them too
		for (int i = 0; i < lib.PrecursorCount; i++)
		{
			MslLibraryEntry skeleton = lib.GetEntry(i, includeFragments: false)!;
			MslLibraryEntry e = expected.Single(x => x.FullSequence == skeleton.FullSequence);
			AssertSameFlags(skeleton, e, "skeleton " + e.FullSequence);
		}
	}

	/// <summary>The bug this change fixes: a read and rewrite erased the bits.</summary>
	[TestCase(false)]
	[TestCase(true)]
	public void ReadAndWriteAgain_KeepsStatusAndStrike_ByteForByte(bool indexOnly)
	{
		string first = TempPath($"first_{indexOnly}");
		string second = TempPath($"second_{indexOnly}");
		MslWriter.Write(first, EveryCombination());

		using (MslLibrary lib = indexOnly ? MslLibrary.LoadIndexOnly(first) : MslLibrary.Load(first))
			MslLibrary.Save(second, lib.GetAllEntries().ToList());

		Assert.That(File.ReadAllBytes(second), Is.EqualTo(File.ReadAllBytes(first)));
	}

	[TestCase(false)]
	[TestCase(true)]
	public void AFileWrittenBy1_0_597_ReadsAsNone(bool indexOnly)
	{
		string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "FileReadingTests",
			"SpectralLibraryTests", "SpectralLibraryData", "msl_written_by_1.0.597.msl");
		using MslLibrary lib = indexOnly ? MslLibrary.LoadIndexOnly(path) : MslLibrary.Load(path);

		List<MslLibraryEntry> entries = lib.GetAllEntries().ToList();
		Assert.That(entries, Is.Not.Empty);
		Assert.That(entries.Select(e => e.RtStatus), Is.All.EqualTo(MslFormat.RetentionTimeStatus.None));
		Assert.That(entries.Select(e => e.RtStruckByPrediction), Is.All.False);
	}

	[Test]
	public void Writers_RefuseAnUndefinedStatus_AndWriteNothing()
	{
		foreach (string writer in new[] { "Write", "WriteStreaming" })
		{
			string path = TempPath("undefined_" + writer);
			List<MslLibraryEntry> entries = EveryCombination();
			entries[3].RtStatus = (MslFormat.RetentionTimeStatus)7;

			Assert.That(() => WriteWith(writer, path, entries, 0), Throws.ArgumentException, writer);
			Assert.That(File.Exists(path), Is.False, writer);
			Assert.That(MslWriter.ValidateEntries(entries), Has.Some.Contains("RtStatus"), writer);
		}
	}

	// ── Merge and UpdateAndSave ───────────────────────────────────────────────

	[TestCase(MslMergeConflictPolicy.KeepFirst)]
	[TestCase(MslMergeConflictPolicy.KeepLast)]
	[TestCase(MslMergeConflictPolicy.KeepLowestQValue)]
	public void Merge_CarriesEachEntrysStatus(MslMergeConflictPolicy policy)
	{
		string a = TempPath($"merge_a_{policy}");
		string b = TempPath($"merge_b_{policy}");
		string merged = TempPath($"merge_out_{policy}");

		// One precursor in both sources, with different statuses; one in each source alone
		MslLibraryEntry sharedA = Entry("SHAREDK", 2, 500.0, 20.0, MslFormat.RetentionTimeStatus.Confirmed);
		MslLibraryEntry sharedB = Entry("SHAREDK", 2, 500.0, 20.0, MslFormat.RetentionTimeStatus.Challenged, struck: true);
		sharedA.QValue = 0.001f;
		sharedB.QValue = 0.0001f;
		MslWriter.Write(a, new List<MslLibraryEntry> { sharedA, Entry("ONLYAK", 2, 450.0, 15.0, MslFormat.RetentionTimeStatus.ObservedOnce) });
		MslWriter.Write(b, new List<MslLibraryEntry> { sharedB, Entry("ONLYBK", 2, 550.0, 25.0, MslFormat.RetentionTimeStatus.Challenged, struck: true) });

		MslMerger.Merge(new[] { a, b }, merged, policy);

		using MslLibrary lib = MslLibrary.Load(merged);
		Assert.That(lib.TryGetEntry("ONLYAK", 2, out MslLibraryEntry? onlyA), Is.True);
		Assert.That(onlyA!.RtStatus, Is.EqualTo(MslFormat.RetentionTimeStatus.ObservedOnce));
		Assert.That(lib.TryGetEntry("ONLYBK", 2, out MslLibraryEntry? onlyB), Is.True);
		Assert.That(onlyB!.RtStatus, Is.EqualTo(MslFormat.RetentionTimeStatus.Challenged));
		Assert.That(onlyB.RtStruckByPrediction, Is.True);

		MslLibraryEntry winner = policy == MslMergeConflictPolicy.KeepFirst ? sharedA : sharedB;
		Assert.That(lib.TryGetEntry("SHAREDK", 2, out MslLibraryEntry? shared), Is.True);
		AssertSameFlags(shared!, winner, $"shared under {policy}");
	}

	private static LibrarySpectrum Spectrum(string seq, int charge, double mz, double? rt, int ionCount)
	{
		var ions = new List<MatchedFragmentIon>();
		for (int i = 1; i <= ionCount; i++)
		{
			var product = new Product(ProductType.y, FragmentationTerminus.C, neutralMass: 0.0,
				fragmentNumber: i, residuePosition: i, neutralLoss: 0.0);
			ions.Add(new MatchedFragmentIon(product, experMz: 100.0 + 10 * i, experIntensity: 1.0, charge: 1));
		}
		return new LibrarySpectrum(seq, mz, charge, ions, rt);
	}

	/// <summary>
	/// UpdateAndSave checks no evidence: an entry whose time it writes gets None; an entry whose
	/// time it keeps keeps its status and strike.
	/// </summary>
	[Test]
	public void UpdateAndSave_KeepsTheStatusOnlyWhereItKeepsTheTime()
	{
		string source = TempPath("update_source");
		string output = TempPath("update_output");
		MslWriter.Write(source, new List<MslLibraryEntry>
		{
			Entry("NEWTIMEK", 2, 400.0, 10.0, MslFormat.RetentionTimeStatus.Confirmed),
			Entry("NOTIMEK", 2, 450.0, 20.0, MslFormat.RetentionTimeStatus.Challenged, struck: true),
			Entry("KEPTK", 2, 500.0, 30.0, MslFormat.RetentionTimeStatus.ObservedOnce),
		});

		using (MslLibrary lib = MslLibrary.Load(source))
		{
			var incoming = new List<LibrarySpectrum>
			{
				Spectrum("NEWTIMEK", 2, 400.0, rt: 12.5, ionCount: 5), // more ions and a time: replaced
				Spectrum("NOTIMEK", 2, 450.0, rt: null, ionCount: 5),  // more ions, no time: time kept
				Spectrum("KEPTK", 2, 500.0, rt: 31.0, ionCount: 1),    // fewer ions: entry retained
				Spectrum("NOVELK", 2, 550.0, rt: 40.0, ionCount: 3),   // not in the library: added
			};
			MslUpdateResult result = lib.UpdateAndSave(incoming, output);
			Assert.That((result.ReplacedCount, result.RetainedCount, result.AddedCount), Is.EqualTo((2, 1, 1)));
		}

		using MslLibrary updated = MslLibrary.Load(output);
		MslLibraryEntry Get(string seq) => updated.TryGetEntry(seq, 2, out MslLibraryEntry? e) ? e! : throw new AssertionException(seq);

		Assert.That(Get("NEWTIMEK").RtStatus, Is.EqualTo(MslFormat.RetentionTimeStatus.None));
		Assert.That(Get("NOTIMEK").RtStatus, Is.EqualTo(MslFormat.RetentionTimeStatus.Challenged));
		Assert.That(Get("NOTIMEK").RtStruckByPrediction, Is.True);
		Assert.That(Get("NOTIMEK").RetentionTime, Is.EqualTo(20.0).Within(1e-4));
		Assert.That(Get("KEPTK").RtStatus, Is.EqualTo(MslFormat.RetentionTimeStatus.ObservedOnce));
		Assert.That(Get("NOVELK").RtStatus, Is.EqualTo(MslFormat.RetentionTimeStatus.None));
	}

	[Test]
	public void LibrarySpectrumConversion_HasNoStatus()
	{
		MslLibraryEntry entry = Entry("PEPTIDEK", 2, 450.0, 25.0, MslFormat.RetentionTimeStatus.Confirmed, struck: true);
		MslLibraryEntry back = MslLibraryEntry.FromLibrarySpectrum(entry.ToLibrarySpectrum())!;

		Assert.That(back.RtStatus, Is.EqualTo(MslFormat.RetentionTimeStatus.None));
		Assert.That(back.RtStruckByPrediction, Is.False);
	}
}
