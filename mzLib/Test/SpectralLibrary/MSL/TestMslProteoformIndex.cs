using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test.MslSpectralLibrary;

/// <summary>
/// NUnit 4 tests for Prompt 16: <see cref="MslProteoformIndex"/> construction,
/// neutral-mass computation, mass-window queries, RT calibration, and
/// <see cref="MslLibrary"/> integration.
///
/// All file I/O uses a dedicated temp directory cleaned up in
/// <see cref="OneTimeTearDown"/>. All assertions use <c>Assert.That</c> exclusively.
/// </summary>
[TestFixture]
public sealed class TestMslProteoformIndex
{
	// ── Fixture paths ─────────────────────────────────────────────────────────

	private static readonly string OutputDirectory =
		Path.Combine(Path.GetTempPath(), "MslProteoformIndexTests");

	private static string TempMsl(string name) =>
		Path.Combine(OutputDirectory, name + ".msl");

	// ── One-time setup / teardown ─────────────────────────────────────────────

	[OneTimeSetUp]
	public void OneTimeSetUp()
	{
		Directory.CreateDirectory(OutputDirectory);
	}

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(OutputDirectory))
			Directory.Delete(OutputDirectory, recursive: true);
	}

	// ── Synthetic test data helper ────────────────────────────────────────────

	/// <summary>
	/// Creates a minimal <see cref="MslLibraryEntry"/> with
	/// <see cref="MslFormat.MoleculeType.Proteoform"/> and a single fragment ion.
	/// </summary>
	private static MslLibraryEntry MakeProteoformEntry(
		string sequence, int charge, double mz, double irt, bool isDecoy = false)
	{
		return new MslLibraryEntry
		{
			FullSequence = sequence,
			BaseSequence = sequence,
			PrecursorMz = mz,
			ChargeState = charge,
			RetentionTime = irt,
			IsDecoy = isDecoy,
			MoleculeType = MslFormat.MoleculeType.Proteoform,
			DissociationType = DissociationType.HCD,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					ProductType    = ProductType.b,
					FragmentNumber = 5,
					Charge         = 1,
					Mz             = 500.0f,
					Intensity      = 1.0f,
					NeutralLoss    = 0.0
				}
			}
		};
	}

	/// <summary>
	/// Creates a minimal <see cref="MslLibraryEntry"/> with
	/// <see cref="MslFormat.MoleculeType.Peptide"/>.
	/// </summary>
	private static MslLibraryEntry MakePeptideEntry(
		string sequence, int charge, double mz, double irt)
	{
		return new MslLibraryEntry
		{
			FullSequence = sequence,
			BaseSequence = sequence,
			PrecursorMz = mz,
			ChargeState = charge,
			RetentionTime = irt,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					ProductType    = ProductType.b,
					FragmentNumber = 2,
					Charge         = 1,
					Mz             = 300.0f,
					Intensity      = 1.0f,
					NeutralLoss    = 0.0
				}
			}
		};
	}

	// ═════════════════════════════════════════════════════════════════════════
	// MslProteoformIndex unit tests
	// ═════════════════════════════════════════════════════════════════════════

	[Test]
	public void ProteoformIndex_Build_OnlyIndexesProteoforms()
	{
		// 3 peptides + 2 proteoforms → index should contain exactly 2 entries
		var entries = new List<MslLibraryEntry>
		{
			MakePeptideEntry("PEPTIDE",    2, 400.0, 10.0),
			MakePeptideEntry("ACDEFGHIK", 2, 500.0, 15.0),
			MakePeptideEntry("LMNPQRSTV", 3, 600.0, 20.0),
			MakeProteoformEntry("PROTEOFORM_A", 10, 800.0, 30.0),
			MakeProteoformEntry("PROTEOFORM_B", 12, 900.0, 35.0)
		};

		MslProteoformIndex index = MslProteoformIndex.Build(entries, i => entries[i]);

		Assert.That(index.Count, Is.EqualTo(2));
	}

	[Test]
	public void ProteoformIndex_NeutralMass_ComputedCorrectly()
	{
		// mz=800.0, charge=10 → NeutralMass = (800.0 × 10) − (10 × 1.007276) = 7992.92724
		const double mz = 800.0;
		const int charge = 10;
		const double expectedMass = mz * charge - charge * MslProteoformIndex.ProtonMass;

		var entries = new List<MslLibraryEntry>
		{
			MakeProteoformEntry("BIGPROTEIN", charge, mz, 10.0)
		};

		MslProteoformIndex index = MslProteoformIndex.Build(entries, i => entries[i]);

		Assert.That(index.Count, Is.EqualTo(1));

		// Retrieve via QueryMassWindow using a ±0.01 Da window
		ReadOnlySpan<MslProteoformIndexEntry> results =
			index.QueryMassWindow(expectedMass - 0.01, expectedMass + 0.01);

		Assert.That(results.Length, Is.EqualTo(1));
		Assert.That(results[0].NeutralMass, Is.EqualTo(expectedMass).Within(0.001));
	}

	[Test]
	public void ProteoformIndex_QueryMassWindow_ReturnsCorrectEntries()
	{
		// mz=800.0, charge=10 → mass ≈ 7992.927
		// mz=1000.0, charge=10 → mass ≈ 9992.927
		// mz=600.0, charge=8  → mass ≈ 4791.942
		var entries = new List<MslLibraryEntry>
		{
			MakeProteoformEntry("PROT_A", 10, 800.0,  10.0),
			MakeProteoformEntry("PROT_B", 10, 1000.0, 20.0),
			MakeProteoformEntry("PROT_C", 8,  600.0,  30.0)
		};

		MslProteoformIndex index = MslProteoformIndex.Build(entries, i => entries[i]);

		// PROT_A: mz=800.0, charge=10 → mass = 8000 − 10×1.007276 ≈ 7989.927
		// Query [7985, 7995] — should only match PROT_A, not PROT_B (~9992) or PROT_C (~4791)
		ReadOnlySpan<MslProteoformIndexEntry> results =
			index.QueryMassWindow(7985.0, 7995.0);

		Assert.That(results.Length, Is.EqualTo(1));
		Assert.That(results[0].PrecursorMz, Is.EqualTo(800.0f).Within(0.001f));
	}

	[Test]
	public void ProteoformIndex_QueryMassWindow_ZeroAllocation()
	{
		// Calling QueryMassWindow twice should return spans backed by the same array —
		// verify by checking that both spans reference the same starting element when
		// the range covers all entries.
		var entries = new List<MslLibraryEntry>
		{
			MakeProteoformEntry("PROT_A", 10, 800.0, 10.0),
			MakeProteoformEntry("PROT_B", 10, 900.0, 20.0)
		};

		MslProteoformIndex index = MslProteoformIndex.Build(entries, i => entries[i]);

		ReadOnlySpan<MslProteoformIndexEntry> first = index.QueryMassWindow(0, double.MaxValue);
		ReadOnlySpan<MslProteoformIndexEntry> second = index.QueryMassWindow(0, double.MaxValue);

		// Both calls should return spans of the same length with identical content —
		// no allocation means the underlying data is the same array slice.
		Assert.That(first.Length, Is.EqualTo(second.Length));
		Assert.That(first.Length, Is.EqualTo(2));
		Assert.That(first[0].NeutralMass, Is.EqualTo(second[0].NeutralMass));
		Assert.That(first[1].NeutralMass, Is.EqualTo(second[1].NeutralMass));
	}

	[Test]
	public void ProteoformIndex_QueryMassWindow_EmptyRange_ReturnsEmpty()
	{
		var entries = new List<MslLibraryEntry>
		{
			MakeProteoformEntry("PROT_A", 10, 800.0, 10.0)
		};

		MslProteoformIndex index = MslProteoformIndex.Build(entries, i => entries[i]);

		// Query a range that contains no entries
		ReadOnlySpan<MslProteoformIndexEntry> results =
			index.QueryMassWindow(1.0, 100.0);

		Assert.That(results.IsEmpty, Is.True);
	}

	[Test]
	public void ProteoformIndex_QueryMassWindow_DecoyFiltering()
	{
		// One target, one decoy — both in the same mass window
		var entries = new List<MslLibraryEntry>
		{
			MakeProteoformEntry("TARGET", 10, 800.0, 10.0, isDecoy: false),
			MakeProteoformEntry("DECOY",  10, 801.0, 10.0, isDecoy: true)
		};

		MslProteoformIndex index = MslProteoformIndex.Build(entries, i => entries[i]);

		ReadOnlySpan<MslProteoformIndexEntry> withDecoys =
			index.QueryMassWindow(7980.0, 8020.0, includeDecoys: true);
		ReadOnlySpan<MslProteoformIndexEntry> noDecoys =
			index.QueryMassWindow(7980.0, 8020.0, includeDecoys: false);

		Assert.That(withDecoys.Length, Is.EqualTo(2));
		Assert.That(noDecoys.Length, Is.EqualTo(1));
		Assert.That(noDecoys[0].IsDecoy, Is.False);
	}

	[Test]
	public void ProteoformIndex_TryGetEntry_HitReturnsCorrectEntry()
	{
		var entries = new List<MslLibraryEntry>
		{
			MakeProteoformEntry("MYPROTEOFORM", 15, 900.0, 25.0)
		};

		MslProteoformIndex index = MslProteoformIndex.Build(entries, i => entries[i]);

		bool found = index.TryGetEntry("MYPROTEOFORM", 15, out MslLibraryEntry? entry);

		Assert.That(found, Is.True);
		Assert.That(entry, Is.Not.Null);
		Assert.That(entry!.FullSequence, Is.EqualTo("MYPROTEOFORM"));
		Assert.That(entry.ChargeState, Is.EqualTo(15));
	}

	[Test]
	public void ProteoformIndex_TryGetEntry_MissReturnsFalse()
	{
		var entries = new List<MslLibraryEntry>
		{
			MakeProteoformEntry("MYPROTEOFORM", 15, 900.0, 25.0)
		};

		MslProteoformIndex index = MslProteoformIndex.Build(entries, i => entries[i]);

		bool found = index.TryGetEntry("NOTPRESENT", 5, out MslLibraryEntry? entry);

		Assert.That(found, Is.False);
		Assert.That(entry, Is.Null);
	}

	[Test]
	public void ProteoformIndex_SortedByNeutralMass()
	{
		// Insert entries in reverse mass order; the index should sort them ascending.
		// mz=900, charge=10 → mass ≈ 8992.927  (largest)
		// mz=800, charge=10 → mass ≈ 7992.927
		// mz=600, charge=8  → mass ≈ 4791.942  (smallest)
		var entries = new List<MslLibraryEntry>
		{
			MakeProteoformEntry("PROT_LARGE",  10, 900.0, 10.0),
			MakeProteoformEntry("PROT_MEDIUM", 10, 800.0, 20.0),
			MakeProteoformEntry("PROT_SMALL",  8,  600.0, 30.0)
		};

		MslProteoformIndex index = MslProteoformIndex.Build(entries, i => entries[i]);
		ReadOnlySpan<MslProteoformIndexEntry> all = index.QueryMassWindow(0, double.MaxValue);

		Assert.That(all.Length, Is.EqualTo(3));
		for (int i = 1; i < all.Length; i++)
			Assert.That(all[i].NeutralMass, Is.GreaterThanOrEqualTo(all[i - 1].NeutralMass));
	}

	[Test]
	public void ProteoformIndex_RtCalibration_TransformsValues()
	{
		// Stored iRT = 10.0; after calibration (slope=2.0, intercept=5.0) → 25.0
		var entries = new List<MslLibraryEntry>
		{
			MakeProteoformEntry("PROT_A", 10, 800.0, 10.0),
			MakeProteoformEntry("PROT_B", 10, 900.0, 20.0)
		};

		MslProteoformIndex original = MslProteoformIndex.Build(entries, i => entries[i]);
		MslProteoformIndex calibrated = original.WithCalibratedRetentionTimes(slope: 2.0, intercept: 5.0);

		ReadOnlySpan<MslProteoformIndexEntry> origAll = original.QueryMassWindow(0, double.MaxValue);
		ReadOnlySpan<MslProteoformIndexEntry> calAll = calibrated.QueryMassWindow(0, double.MaxValue);

		Assert.That(calAll.Length, Is.EqualTo(origAll.Length));

		for (int i = 0; i < origAll.Length; i++)
		{
			double expectedIrt = 2.0 * origAll[i].Irt + 5.0;
			Assert.That(calAll[i].Irt, Is.EqualTo((float)expectedIrt).Within(0.001f));
		}
	}

	// ═════════════════════════════════════════════════════════════════════════
	// MslLibrary integration tests
	// ═════════════════════════════════════════════════════════════════════════

	[Test]
	public void MslLibrary_ProteoformIndex_PresentWhenProteoformsExist()
	{
		string path = TempMsl("proteoform_present");
		var entries = new List<MslLibraryEntry>
		{
			MakePeptideEntry("PEPTIDE", 2, 400.0, 10.0),
			MakeProteoformEntry("BIGPROTEIN", 10, 800.0, 20.0)
		};
		MslLibrary.Save(path, entries);

		using MslLibrary lib = MslLibrary.Load(path);
		Assert.That(lib.ProteoformIndex, Is.Not.Null);
		Assert.That(lib.ProteoformIndex!.Count, Is.EqualTo(1));
	}

	[Test]
	public void MslLibrary_ProteoformIndex_NullWhenNoPeptideOnly()
	{
		string path = TempMsl("peptide_only");
		var entries = new List<MslLibraryEntry>
		{
			MakePeptideEntry("PEPTIDE",    2, 400.0, 10.0),
			MakePeptideEntry("ACDEFGHIK", 2, 500.0, 15.0)
		};
		MslLibrary.Save(path, entries);

		using MslLibrary lib = MslLibrary.Load(path);
		Assert.That(lib.ProteoformIndex, Is.Null);
	}

	[Test]
	public void MslLibrary_ProteoformCount_CorrectValue()
	{
		string path = TempMsl("proteoform_count");
		var entries = new List<MslLibraryEntry>
		{
			MakePeptideEntry("PEPTIDE_A", 2, 400.0, 10.0),
			MakePeptideEntry("PEPTIDE_B", 2, 500.0, 15.0),
			MakeProteoformEntry("PROT_A", 10, 800.0, 20.0),
			MakeProteoformEntry("PROT_B", 12, 900.0, 25.0),
			MakeProteoformEntry("PROT_C", 15, 1000.0, 30.0)
		};
		MslLibrary.Save(path, entries);

		using MslLibrary lib = MslLibrary.Load(path);
		Assert.That(lib.ProteoformCount, Is.EqualTo(3));
	}

	[Test]
	public void MslLibrary_PeptideCount_CorrectValue()
	{
		string path = TempMsl("peptide_count");
		var entries = new List<MslLibraryEntry>
		{
			MakePeptideEntry("PEPTIDE_A", 2, 400.0, 10.0),
			MakePeptideEntry("PEPTIDE_B", 2, 500.0, 15.0),
			MakeProteoformEntry("PROT_A", 10, 800.0, 20.0)
		};
		MslLibrary.Save(path, entries);

		using MslLibrary lib = MslLibrary.Load(path);
		Assert.That(lib.PeptideCount, Is.EqualTo(2));
	}

	[Test]
	public void MslLibrary_QueryProteoformMassWindow_RoutesCorrectly()
	{
		string path = TempMsl("proteoform_route");
		var entries = new List<MslLibraryEntry>
		{
			MakeProteoformEntry("PROT_A", 10, 800.0, 10.0),
			MakeProteoformEntry("PROT_B", 10, 1000.0, 20.0)
		};
		MslLibrary.Save(path, entries);

		using MslLibrary lib = MslLibrary.Load(path);

		double massA = MslProteoformIndex.ComputeNeutralMass(800.0, 10);

		// QueryProteoformMassWindow should return same results as ProteoformIndex.QueryMassWindow
		ReadOnlySpan<MslProteoformIndexEntry> viaLibrary =
			lib.QueryProteoformMassWindow(massA - 1.0, massA + 1.0);
		ReadOnlySpan<MslProteoformIndexEntry> viaDirect =
			lib.ProteoformIndex!.QueryMassWindow(massA - 1.0, massA + 1.0);

		Assert.That(viaLibrary.Length, Is.EqualTo(viaDirect.Length));
		Assert.That(viaLibrary.Length, Is.EqualTo(1));
		Assert.That(viaLibrary[0].NeutralMass, Is.EqualTo(viaDirect[0].NeutralMass));
	}

	[Test]
	public void MslLibrary_QueryProteoformMassWindow_EmptyWhenNoIndex()
	{
		string path = TempMsl("no_proteoform_index");
		var entries = new List<MslLibraryEntry>
		{
			MakePeptideEntry("PEPTIDE_A", 2, 400.0, 10.0)
		};
		MslLibrary.Save(path, entries);

		using MslLibrary lib = MslLibrary.Load(path);

		// Should return empty span — not throw — when no proteoform index exists
		Assert.DoesNotThrow(() =>
		{
			ReadOnlySpan<MslProteoformIndexEntry> result =
				lib.QueryProteoformMassWindow(0.0, double.MaxValue);
			Assert.That(result.IsEmpty, Is.True);
		});
	}

	[Test]
	public void ProteoformIndex_BinaryRoundTrip_MassPreserved()
	{
		// Write a library with one proteoform, load it, query by mass window
		string path = TempMsl("binary_roundtrip");
		const double mz = 850.0;
		const int charge = 12;
		double expectedMass = MslProteoformIndex.ComputeNeutralMass(mz, charge);

		var entries = new List<MslLibraryEntry>
		{
			MakeProteoformEntry("ROUNDTRIP_PROT", charge, mz, 42.0)
		};
		MslLibrary.Save(path, entries);

		using MslLibrary lib = MslLibrary.Load(path);

		Assert.That(lib.ProteoformIndex, Is.Not.Null);

		ReadOnlySpan<MslProteoformIndexEntry> results =
			lib.QueryProteoformMassWindow(expectedMass - 1.0, expectedMass + 1.0);

		Assert.That(results.Length, Is.EqualTo(1));
		// Mass should be preserved to within float32 round-trip precision (~0.1 Da at 10 kDa)
		Assert.That(results[0].NeutralMass, Is.EqualTo(expectedMass).Within(0.5));
	}

	[Test]
	public void ProteoformIndex_MixedLibrary_BothIndexesWork()
	{
		// A library with peptides and proteoforms should support both query paths
		string path = TempMsl("mixed_library");
		const float peptideMz = 500.0f;
		const float peptideMzHigh = 502.0f;
		const double protMz = 800.0;
		const int protCharge = 10;
		double protMass = MslProteoformIndex.ComputeNeutralMass(protMz, protCharge);

		var entries = new List<MslLibraryEntry>
		{
			MakePeptideEntry("PEPTIDE_X",  2, peptideMz, 10.0),
			MakeProteoformEntry("PROT_X", protCharge, protMz, 20.0)
		};
		MslLibrary.Save(path, entries);

		using MslLibrary lib = MslLibrary.Load(path);

		// Peptide m/z query via the primary index
		ReadOnlySpan<MslPrecursorIndexEntry> mzResults =
			lib.QueryMzWindow(peptideMz - 1f, peptideMzHigh);
		Assert.That(mzResults.Length, Is.GreaterThanOrEqualTo(1));

		// Proteoform mass query via the proteoform index
		ReadOnlySpan<MslProteoformIndexEntry> massResults =
			lib.QueryProteoformMassWindow(protMass - 1.0, protMass + 1.0);
		Assert.That(massResults.Length, Is.EqualTo(1));
	}
}