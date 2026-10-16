using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Omics.SpectrumMatch;
using Readers.SpectralLibrary;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test.MslSpectralLibrary;

/// <summary>
/// Tests for the precursor-level entrapment flag (PrecursorFlags bit 3) and the
/// rt_is_calibrated flag (bit 2), plus a modification-position round-trip.
///
/// Entrapment is independent of decoy status, so the four combinations T / D / ET / ED must
/// each survive write → read, index-only load, merge and RT calibration. The bit was added
/// without a format-version bump because readers that predate it ignore unknown bits.
/// </summary>
[TestFixture]
public sealed class TestMslEntrapmentAndRtCalibrated
{
	private readonly List<string> _tempPaths = new();

	private string NewTempMsl()
	{
		string raw = Path.GetTempFileName();
		string path = Path.ChangeExtension(raw, ".msl");
		File.Delete(raw);
		_tempPaths.Add(path);
		return path;
	}

	[TearDown]
	public void TearDown()
	{
		foreach (string p in _tempPaths)
			if (File.Exists(p)) File.Delete(p);
		_tempPaths.Clear();
	}

	private static MslLibraryEntry MakeEntry(string seq, double mz, double rt,
		bool isDecoy, bool isEntrapment, bool rtIsCalibrated = false, string? fullSeq = null)
	{
		return new MslLibraryEntry
		{
			FullSequence = fullSeq ?? seq,
			BaseSequence = seq,
			PrecursorMz = mz,
			ChargeState = 2,
			RetentionTime = rt,
			RtIsCalibrated = rtIsCalibrated,
			IsDecoy = isDecoy,
			IsEntrapment = isEntrapment,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			QValue = float.NaN,
			ProteinAccession = string.Empty,
			ProteinName = string.Empty,
			GeneName = string.Empty,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz = 300.0f, Intensity = 1.0f, ProductType = ProductType.y,
					FragmentNumber = 2, ResiduePosition = 2, Charge = 1, NeutralLoss = 0.0
				}
			}
		};
	}

	/// <summary>Target, decoy, entrapment target and entrapment decoy, in that order.</summary>
	private static MslLibraryEntry[] FourLabels() => new[]
	{
		MakeEntry("TARGETK", 500.0, 10.0, isDecoy: false, isEntrapment: false),
		MakeEntry("DECOYK", 501.0, 20.0, isDecoy: true, isEntrapment: false),
		MakeEntry("ENTRAPTK", 502.0, 30.0, isDecoy: false, isEntrapment: true),
		MakeEntry("ENTRAPDK", 503.0, 40.0, isDecoy: true, isEntrapment: true, rtIsCalibrated: true),
	};

	private static void AssertFourLabels(IEnumerable<MslLibraryEntry> loaded)
	{
		var bySeq = loaded.ToDictionary(e => e.BaseSequence);
		Assert.That(bySeq, Has.Count.EqualTo(4));
		Assert.That((bySeq["TARGETK"].IsDecoy, bySeq["TARGETK"].IsEntrapment), Is.EqualTo((false, false)));
		Assert.That((bySeq["DECOYK"].IsDecoy, bySeq["DECOYK"].IsEntrapment), Is.EqualTo((true, false)));
		Assert.That((bySeq["ENTRAPTK"].IsDecoy, bySeq["ENTRAPTK"].IsEntrapment), Is.EqualTo((false, true)));
		Assert.That((bySeq["ENTRAPDK"].IsDecoy, bySeq["ENTRAPDK"].IsEntrapment), Is.EqualTo((true, true)));
		Assert.That(bySeq["ENTRAPDK"].RtIsCalibrated, Is.True);
		Assert.That(bySeq.Values.Count(e => e.RtIsCalibrated), Is.EqualTo(1));
	}

	// ── Flag encoding ─────────────────────────────────────────────────────────

	[Test]
	public void EncodeDecode_AllSixteenCombinations_RoundTrip()
	{
		foreach (bool d in new[] { false, true })
		foreach (bool p in new[] { false, true })
		foreach (bool r in new[] { false, true })
		foreach (bool e in new[] { false, true })
		{
			byte flags = MslFormat.EncodePrecursorFlags(d, p, r, e);
			Assert.That(MslFormat.DecodePrecursorFlags(flags), Is.EqualTo((d, p, r)));
			Assert.That(MslFormat.DecodeIsEntrapment(flags), Is.EqualTo(e));
			Assert.That(flags & 0b1111_0000, Is.Zero, "bits 4-7 must stay reserved");
		}
	}

	[Test]
	public void Encode_EntrapmentIsBit3_AndThreeArgOverloadLeavesItClear()
	{
		Assert.That(MslFormat.EncodePrecursorFlags(false, false, false, isEntrapment: true), Is.EqualTo(0b_0000_1000));
		Assert.That(MslFormat.EncodePrecursorFlags(true, true, true),
			Is.EqualTo(MslFormat.EncodePrecursorFlags(true, true, true, isEntrapment: false)));
	}

	[Test]
	public void Decode_EntrapmentBitDoesNotChangeDecoyReading()
	{
		// A reader that predates bit 3 sees an entrapment target as a target and an
		// entrapment decoy as a decoy — the documented backward-compatible fallback.
		byte et = MslFormat.EncodePrecursorFlags(false, false, false, isEntrapment: true);
		byte ed = MslFormat.EncodePrecursorFlags(true, false, false, isEntrapment: true);
		Assert.That(MslFormat.DecodePrecursorFlags(et).isDecoy, Is.False);
		Assert.That(MslFormat.DecodePrecursorFlags(ed).isDecoy, Is.True);
	}

	// ── Round trips ───────────────────────────────────────────────────────────

	[Test]
	public void SaveLoad_PreservesEntrapmentAndRtCalibrated()
	{
		string path = NewTempMsl();
		MslLibrary.Save(path, FourLabels());

		using MslLibrary lib = MslLibrary.Load(path);
		AssertFourLabels(lib.GetAllEntries(includeDecoys: true));
	}

	[Test]
	public void WriteStreaming_PreservesEntrapmentAndRtCalibrated()
	{
		string path = NewTempMsl();
		MslWriter.WriteStreaming(path, FourLabels());

		using MslLibrary lib = MslLibrary.Load(path);
		AssertFourLabels(lib.GetAllEntries(includeDecoys: true));
	}

	[Test]
	public void LoadIndexOnly_PreservesEntrapmentAndRtCalibrated()
	{
		string path = NewTempMsl();
		MslLibrary.Save(path, FourLabels());

		using MslLibrary lib = MslLibrary.LoadIndexOnly(path);
		AssertFourLabels(lib.GetAllEntries(includeDecoys: true));
	}

	[Test]
	public void Merge_PreservesEntrapmentAndRtCalibrated()
	{
		MslLibraryEntry[] all = FourLabels();
		string a = NewTempMsl(), b = NewTempMsl(), merged = NewTempMsl();
		MslLibrary.Save(a, all.Take(2).ToArray());
		MslLibrary.Save(b, all.Skip(2).ToArray());

		MslMerger.Merge(new[] { a, b }, merged);

		using MslLibrary lib = MslLibrary.Load(merged);
		AssertFourLabels(lib.GetAllEntries(includeDecoys: true));
	}

	// ── Index and QueryWindow ─────────────────────────────────────────────────

	[Test]
	public void QueryWindow_EntrapmentTargetsCompeteAsTargets_EntrapmentDecoysFollowDecoyFilter()
	{
		string path = NewTempMsl();
		MslLibrary.Save(path, FourLabels());
		using MslLibrary lib = MslLibrary.LoadIndexOnly(path);

		using (MslWindowResults targetsOnly = lib.QueryWindow(400f, 600f, 0f, 100f, includeDecoys: false))
		{
			var entries = targetsOnly.Entries.ToArray();
			Assert.That(entries.Length, Is.EqualTo(2));
			Assert.That(entries.Count(e => e.IsEntrapment), Is.EqualTo(1), "ET is returned with the targets");
			Assert.That(entries.All(e => e.IsDecoy == 0));
		}

		using (MslWindowResults all = lib.QueryWindow(400f, 600f, 0f, 100f, includeDecoys: true))
		{
			var entries = all.Entries.ToArray();
			Assert.That(entries.Length, Is.EqualTo(4));
			Assert.That(entries.Count(e => e.IsEntrapment), Is.EqualTo(2));
			Assert.That(entries.Count(e => e.IsEntrapment && e.IsDecoy != 0), Is.EqualTo(1));
			Assert.That(entries.All(e => e.Flags >> 3 == 0), "index flag bits 3-7 stay reserved");
		}
	}

	[Test]
	public void WithCalibratedRetentionTimes_PreservesEntrapmentFlag()
	{
		string path = NewTempMsl();
		MslLibrary.Save(path, FourLabels());
		using MslLibrary lib = MslLibrary.Load(path);
		using MslLibrary calibrated = lib.WithCalibratedRetentionTimes(slope: 2.0, intercept: 1.0);

		// Original RTs 10/20/30/40 map to 21/41/61/81
		using MslWindowResults hits = calibrated.QueryWindow(400f, 600f, 55f, 90f, includeDecoys: true);
		var entries = hits.Entries.ToArray();
		Assert.That(entries.Length, Is.EqualTo(2));
		Assert.That(entries.All(e => e.IsEntrapment), Is.True);
	}

	// ── UpdateAndSave: rt_is_calibrated reflects whether RT was normalised ────

	private static LibrarySpectrum MakeSpectrum(string seq, double precMz, double observedRt, int fragCount)
	{
		var ions = new List<MatchedFragmentIon>(fragCount);
		for (int i = 1; i <= fragCount; i++)
		{
			var product = new Product(ProductType.y, FragmentationTerminus.C, neutralMass: 0.0,
				fragmentNumber: i, residuePosition: i, neutralLoss: 0.0);
			ions.Add(new MatchedFragmentIon(product, experMz: 100.0 + i * 10, experIntensity: 1.0f, charge: 1));
		}
		return new LibrarySpectrum(seq, precMz, 2, ions, observedRt);
	}

	[Test]
	public void UpdateAndSave_WithoutNormalisation_MarksNewEntriesRtCalibrated_AndKeepsEntrapment()
	{
		string path = NewTempMsl(), outPath = NewTempMsl();
		MslLibrary.Save(path, FourLabels());
		using MslLibrary lib = MslLibrary.Load(path);

		// One overlapping spectrum (< minRegressionAnchors) → raw run minutes are stored.
		var incoming = new[]
		{
			MakeSpectrum("ENTRAPTK", 502.0, observedRt: 12.5, fragCount: 5), // replaces ET (more ions)
			MakeSpectrum("NOVELK", 510.0, observedRt: 14.0, fragCount: 3),
		};
		MslUpdateResult result = lib.UpdateAndSave(incoming, outPath);
		Assert.That(result.RtNormalisationApplied, Is.False);

		using MslLibrary updated = MslLibrary.Load(outPath);
		var bySeq = updated.GetAllEntries(includeDecoys: true).ToDictionary(e => e.BaseSequence);
		Assert.That(bySeq["ENTRAPTK"].IsEntrapment, Is.True, "replacement keeps the entrapment label");
		Assert.That(bySeq["ENTRAPTK"].RtIsCalibrated, Is.True);
		Assert.That(bySeq["NOVELK"].RtIsCalibrated, Is.True);
		Assert.That(bySeq["TARGETK"].RtIsCalibrated, Is.False, "untouched entries keep iRT");
	}

	[Test]
	public void UpdateAndSave_WithNormalisation_KeepsNewEntriesOnIrtScale()
	{
		string path = NewTempMsl(), outPath = NewTempMsl();
		MslLibrary.Save(path, FourLabels());
		using MslLibrary lib = MslLibrary.Load(path);

		var incoming = new[]
		{
			MakeSpectrum("TARGETK", 500.0, observedRt: 5.0, fragCount: 1),
			MakeSpectrum("DECOYK", 501.0, observedRt: 10.0, fragCount: 1),
			MakeSpectrum("NOVELK", 510.0, observedRt: 7.5, fragCount: 3),
		};
		MslUpdateResult result = lib.UpdateAndSave(incoming, outPath, minRegressionAnchors: 2);
		Assert.That(result.RtNormalisationApplied, Is.True);

		using MslLibrary updated = MslLibrary.Load(outPath);
		MslLibraryEntry novel = updated.GetAllEntries(includeDecoys: true).Single(e => e.BaseSequence == "NOVELK");
		Assert.That(novel.RtIsCalibrated, Is.False);
		Assert.That(novel.RetentionTime, Is.EqualTo(15.0).Within(1e-3), "7.5 min maps onto iRT 15");
	}

	[Test]
	public void UpdateAndSave_WithNormalisation_LibraryInMinutes_NewEntriesStayInMinutes()
	{
		// Every library entry holds run minutes (RtIsCalibrated = true), e.g. an empirical library.
		// The regression maps observed RT onto that scale, so the new entries hold minutes too.
		string path = NewTempMsl(), outPath = NewTempMsl();
		MslLibrary.Save(path, FourLabels().Select(e => { e.RtIsCalibrated = true; return e; }).ToArray());
		using MslLibrary lib = MslLibrary.Load(path);

		var incoming = new[]
		{
			MakeSpectrum("TARGETK", 500.0, observedRt: 5.0, fragCount: 3), // replaces TARGETK (more ions)
			MakeSpectrum("DECOYK", 501.0, observedRt: 10.0, fragCount: 1),
			MakeSpectrum("NOVELK", 510.0, observedRt: 7.5, fragCount: 3),
		};
		MslUpdateResult result = lib.UpdateAndSave(incoming, outPath, minRegressionAnchors: 2);
		Assert.That(result.RtNormalisationApplied, Is.True);

		using MslLibrary updated = MslLibrary.Load(outPath);
		var bySeq = updated.GetAllEntries(includeDecoys: true).ToDictionary(e => e.BaseSequence);
		Assert.That(bySeq["TARGETK"].Source, Is.EqualTo(MslFormat.SourceType.Empirical), "TARGETK was replaced");
		Assert.That(bySeq["TARGETK"].RtIsCalibrated, Is.True, "replacement keeps the library's scale");
		Assert.That(bySeq["NOVELK"].RtIsCalibrated, Is.True, "novel entry takes the anchors' scale");
		Assert.That(bySeq["NOVELK"].RetentionTime, Is.EqualTo(15.0).Within(1e-3));
		Assert.That(updated.GetAllEntries(includeDecoys: true).All(e => e.RtIsCalibrated), Is.True);
	}

	[Test]
	public void UpdateAndSave_AnchorsOnDifferentScales_DoesNotNormalise()
	{
		// TARGETK is on iRT, ENTRAPDK in minutes: a line through both is meaningless.
		string path = NewTempMsl(), outPath = NewTempMsl();
		MslLibrary.Save(path, FourLabels());
		using MslLibrary lib = MslLibrary.Load(path);

		var incoming = new[]
		{
			MakeSpectrum("TARGETK", 500.0, observedRt: 5.0, fragCount: 1),
			MakeSpectrum("ENTRAPDK", 503.0, observedRt: 10.0, fragCount: 1),
			MakeSpectrum("NOVELK", 510.0, observedRt: 7.5, fragCount: 3),
		};
		MslUpdateResult result = lib.UpdateAndSave(incoming, outPath, minRegressionAnchors: 2);
		Assert.That(result.AnchorCount, Is.EqualTo(2));
		Assert.That(result.RtNormalisationApplied, Is.False);

		using MslLibrary updated = MslLibrary.Load(outPath);
		MslLibraryEntry novel = updated.GetAllEntries(includeDecoys: true).Single(e => e.BaseSequence == "NOVELK");
		Assert.That(novel.RtIsCalibrated, Is.True);
		Assert.That(novel.RetentionTime, Is.EqualTo(7.5).Within(1e-3), "raw run minutes are stored");
	}

	// ── Modification position survives the round trip ─────────────────────────

	[TestCase("AAAQC[Common Fixed:Carbamidomethyl on C]YIDLIIK", "AAAQCYIDLIIK", 'C', 4)]
	[TestCase("PEPM[Common Variable:Oxidation on M]TIDEK", "PEPMTIDEK", 'M', 3)]
	public void SaveLoad_ModifiedSequence_ModStaysOnItsResidue(string fullSeq, string baseSeq,
		char modifiedResidue, int zeroBasedResidueIndex)
	{
		string path = NewTempMsl();
		MslLibrary.Save(path, new[] { MakeEntry(baseSeq, 600.0, 25.0, false, false, fullSeq: fullSeq) });

		using MslLibrary lib = MslLibrary.Load(path);
		Assert.That(lib.TryGetEntry(fullSeq, 2, out MslLibraryEntry? loaded), Is.True);
		Assert.That(loaded!.FullSequence, Is.EqualTo(fullSeq));
		Assert.That(loaded.BaseSequence, Is.EqualTo(baseSeq));

		// The bracket must follow the modified residue, not precede it.
		int bracket = loaded.FullSequence.IndexOf('[');
		Assert.That(loaded.FullSequence[bracket - 1], Is.EqualTo(modifiedResidue));
		Assert.That(loaded.BaseSequence[zeroBasedResidueIndex], Is.EqualTo(modifiedResidue));
		Assert.That(bracket, Is.EqualTo(zeroBasedResidueIndex + 1));
	}
}
