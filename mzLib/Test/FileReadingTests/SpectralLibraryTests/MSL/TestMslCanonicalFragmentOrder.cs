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
/// The writer puts fragments of equal m/z in a fixed order, so the same content always gives the
/// same bytes: whatever order the fragments arrive in, through either writer, and when a file is
/// read and written again.
/// </summary>
[TestFixture]
public sealed class TestMslCanonicalFragmentOrder
{
	private static readonly string OutputDirectory =
		Path.Combine(Path.GetTempPath(), "MslCanonicalFragmentOrderTests");

	[OneTimeSetUp]
	public void OneTimeSetUp() => Directory.CreateDirectory(OutputDirectory);

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(OutputDirectory))
			Directory.Delete(OutputDirectory, recursive: true);
	}

	private static string TempPath(string name) => Path.Combine(OutputDirectory, name + ".msl");

	private static MslFragmentIon Frag(float mz, float intensity, ProductType type, int number, int charge,
		ProductType? secondaryType = null, int secondaryNumber = 0, double neutralLoss = 0.0,
		int residuePosition = 0, bool excludeFromQuant = false) => new()
	{
		Mz = mz,
		Intensity = intensity,
		ProductType = type,
		FragmentNumber = number,
		Charge = charge,
		SecondaryProductType = secondaryType,
		SecondaryFragmentNumber = secondaryNumber,
		NeutralLoss = neutralLoss,
		ResiduePosition = residuePosition,
		ExcludeFromQuant = excludeFromQuant
	};

	/// <summary>
	/// Fragments with several exact m/z ties. b4 1+ and b8 2+ of a poly-alanine peptide have the
	/// same m/z (the case found in a 5.24 M-precursor library); the other ties differ only in a
	/// later field, so every level of the tie-break is needed to order them.
	/// </summary>
	private static List<MslFragmentIon> TiedFragments() => new()
	{
		Frag(285.15573f, 0.206f, ProductType.b, 4, 1, residuePosition: 4),
		Frag(285.15573f, 0.149f, ProductType.b, 8, 2, residuePosition: 8),
		Frag(400.2f, 0.5f, ProductType.y, 3, 1, residuePosition: 13),
		Frag(400.2f, 0.4f, ProductType.y, 3, 1, ProductType.b, 9, residuePosition: 6),
		Frag(500.3f, 1.0f, ProductType.y, 4, 1, neutralLoss: 18.010565, residuePosition: 12),
		Frag(500.3f, 0.9f, ProductType.y, 4, 1, neutralLoss: 17.026549, residuePosition: 12),
		Frag(600.4f, 0.3f, ProductType.y, 5, 1, residuePosition: 11, excludeFromQuant: true),
		Frag(600.4f, 0.3f, ProductType.y, 5, 1, residuePosition: 11, excludeFromQuant: false),
	};

	private static List<MslLibraryEntry> Entries(Func<List<MslFragmentIon>, List<MslFragmentIon>> order) => new()
	{
		new MslLibraryEntry
		{
			FullSequence = "AAAAAAAAAAAAAAAASAGGK",
			BaseSequence = "AAAAAAAAAAAAAAAASAGGK",
			PrecursorMz = 573.3,
			ChargeState = 3,
			RetentionTime = 40.0,
			ProteinAccession = "P00001",
			Source = MslFormat.SourceType.Predicted,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			Nce = 27,
			MatchedFragmentIons = order(TiedFragments())
		}
	};

	private static List<MslFragmentIon> AsGiven(List<MslFragmentIon> f) => f;
	private static List<MslFragmentIon> Reversed(List<MslFragmentIon> f) { f.Reverse(); return f; }
	// A fixed permutation of the 8 fragments (5 is coprime with 8).
	private static List<MslFragmentIon> Shuffled(List<MslFragmentIon> f) =>
		f.Select((x, i) => (x, k: (i * 5 + 3) % f.Count)).OrderBy(t => t.k).Select(t => t.x).ToList();

	private static void WriteWith(string writer, string path, List<MslLibraryEntry> entries, int compressionLevel)
	{
		if (writer == "Write")
			MslWriter.Write(path, entries, compressionLevel);
		else
			MslWriter.WriteStreaming(path, entries, compressionLevel);
	}

	[TestCase("Write", 0)]
	[TestCase("Write", 3)]
	[TestCase("WriteStreaming", 0)]
	[TestCase("WriteStreaming", 3)]
	public void TiedFragments_InAnyInputOrder_WriteIdenticalFiles(string writer, int compressionLevel)
	{
		string a = TempPath($"{writer}_{compressionLevel}_given");
		string b = TempPath($"{writer}_{compressionLevel}_reversed");
		string c = TempPath($"{writer}_{compressionLevel}_shuffled");

		WriteWith(writer, a, Entries(AsGiven), compressionLevel);
		WriteWith(writer, b, Entries(Reversed), compressionLevel);
		WriteWith(writer, c, Entries(Shuffled), compressionLevel);

		byte[] bytesA = File.ReadAllBytes(a);
		Assert.That(File.ReadAllBytes(b), Is.EqualTo(bytesA));
		Assert.That(File.ReadAllBytes(c), Is.EqualTo(bytesA));
	}

	[TestCase("Write")]
	[TestCase("WriteStreaming")]
	public void ReadAndWriteAgain_GivesTheSameBytes(string writer)
	{
		string first = TempPath($"{writer}_first");
		string second = TempPath($"{writer}_second");
		WriteWith(writer, first, Entries(Reversed), 0);

		List<MslLibraryEntry> reread;
		using (MslLibrary lib = MslLibrary.Load(first))
			reread = lib.GetAllEntries().ToList();
		WriteWith(writer, second, reread, 0);

		Assert.That(File.ReadAllBytes(second), Is.EqualTo(File.ReadAllBytes(first)));
	}

	/// <summary>
	/// List.Sort uses insertion sort (stable) below 17 elements and introsort (unstable) above,
	/// so the reordering of tied fragments seen in real libraries needs a longer fragment list.
	/// Here 40 fragments form 20 pairs of exactly equal m/z.
	/// </summary>
	[TestCase("Write")]
	[TestCase("WriteStreaming")]
	public void ManyTiedFragments_ReadAndWriteAgain_GivesTheSameBytes(string writer)
	{
		var frags = new List<MslFragmentIon>();
		for (int i = 0; i < 20; i++)
		{
			float mz = 200f + 37.5f * i;
			frags.Add(Frag(mz, 0.1f + 0.02f * i, ProductType.b, i + 2, 1, residuePosition: i + 2));
			frags.Add(Frag(mz, 0.9f - 0.02f * i, ProductType.b, 2 * (i + 2), 2, residuePosition: 2 * (i + 2)));
		}
		var entries = Entries(_ => frags);

		string first = TempPath($"{writer}_many_first");
		string second = TempPath($"{writer}_many_second");
		WriteWith(writer, first, entries, 0);

		List<MslLibraryEntry> reread;
		using (MslLibrary lib = MslLibrary.Load(first))
			reread = lib.GetAllEntries().ToList();
		WriteWith(writer, second, reread, 0);

		Assert.That(File.ReadAllBytes(second), Is.EqualTo(File.ReadAllBytes(first)));
	}

	/// <summary>
	/// A custom loss is stored as an index into the file's custom-loss table. The indices are
	/// handed out after the fragments are sorted, so two custom losses tied on m/z get the same
	/// indices whatever order the fragments arrive in.
	/// </summary>
	[TestCase("Write")]
	[TestCase("WriteStreaming")]
	public void TiedCustomLosses_InAnyInputOrder_WriteIdenticalFiles(string writer)
	{
		static List<MslFragmentIon> CustomLossFragments() => new()
		{
			Frag(500.3f, 1.0f, ProductType.y, 4, 1, neutralLoss: 123.4567, residuePosition: 12),
			Frag(500.3f, 0.9f, ProductType.y, 4, 1, neutralLoss: 98.7654, residuePosition: 12),
		};

		string given = TempPath($"{writer}_custom_given");
		string reversed = TempPath($"{writer}_custom_reversed");
		string again = TempPath($"{writer}_custom_again");
		WriteWith(writer, given, Entries(_ => CustomLossFragments()), 0);
		WriteWith(writer, reversed, Entries(_ => Reversed(CustomLossFragments())), 0);

		List<MslLibraryEntry> reread;
		using (MslLibrary lib = MslLibrary.Load(given))
			reread = lib.GetAllEntries().ToList();
		WriteWith(writer, again, reread, 0);

		Assert.That(reread.Single().MatchedFragmentIons.Select(f => f.NeutralLoss),
			Is.EqualTo(new[] { 98.7654, 123.4567 }).Within(1e-4));
		byte[] bytesGiven = File.ReadAllBytes(given);
		Assert.That(File.ReadAllBytes(reversed), Is.EqualTo(bytesGiven));
		Assert.That(File.ReadAllBytes(again), Is.EqualTo(bytesGiven));
	}

	[Test]
	public void WrittenOrder_IsMzThenEachTieBreakField()
	{
		string path = TempPath("order");
		MslWriter.Write(path, Entries(Reversed));

		List<MslFragmentIon> f;
		using (MslLibrary lib = MslLibrary.Load(path))
			f = lib.GetAllEntries().Single().MatchedFragmentIons;

		Assert.That(f.Select(x => x.Mz), Is.Ordered);
		// b4 1+ before b8 2+ (fragment number), terminal y3 before internal y3-b9 (secondary type),
		// smaller neutral-loss value first (neutral loss), quantified before excluded.
		Assert.That(f[0].FragmentNumber, Is.EqualTo(4));
		Assert.That(f[1].FragmentNumber, Is.EqualTo(8));
		Assert.That(f[2].SecondaryProductType, Is.Null);
		Assert.That(f[3].SecondaryProductType, Is.EqualTo(ProductType.b));
		Assert.That(f[4].NeutralLoss, Is.LessThan(f[5].NeutralLoss));
		Assert.That(f[6].ExcludeFromQuant, Is.False);
		Assert.That(f[7].ExcludeFromQuant, Is.True);
	}

	[Test]
	public void Comparer_BreaksTiesOnEachStoredField()
	{
		MslFragmentIon baseline = Frag(300f, 0.5f, ProductType.b, 3, 1, ProductType.y, 5, -18.0, 2);
		Assert.That(MslWriter.CompareFragmentsForWrite(baseline, Frag(300f, 0.5f, ProductType.b, 3, 1, ProductType.y, 5, -18.0, 2)), Is.Zero);

		void Less(MslFragmentIon smaller)
		{
			Assert.That(MslWriter.CompareFragmentsForWrite(smaller, baseline), Is.Negative);
			Assert.That(MslWriter.CompareFragmentsForWrite(baseline, smaller), Is.Positive);
		}

		Less(Frag(299f, 0.5f, ProductType.b, 3, 1, ProductType.y, 5, -18.0, 2));           // m/z
		Less(Frag(300f, 0.5f, ProductType.a, 3, 1, ProductType.y, 5, -18.0, 2));           // ion type
		Less(Frag(300f, 0.5f, ProductType.b, 2, 1, ProductType.y, 5, -18.0, 2));           // number
		Less(Frag(300f, 0.5f, ProductType.b, 3, 0, ProductType.y, 5, -18.0, 2));           // charge
		Less(Frag(300f, 0.5f, ProductType.b, 3, 1, null, 5, -18.0, 2));                    // terminal first
		Less(Frag(300f, 0.5f, ProductType.b, 3, 1, ProductType.y, 4, -18.0, 2));           // end residue
		Less(Frag(300f, 0.5f, ProductType.b, 3, 1, ProductType.y, 5, -19.0, 2));           // neutral loss
		Less(Frag(300f, 0.5f, ProductType.b, 3, 1, ProductType.y, 5, -18.0, 1));           // residue
		Less(Frag(300f, 0.4f, ProductType.b, 3, 1, ProductType.y, 5, -18.0, 2));           // intensity
		Assert.That(MslWriter.CompareFragmentsForWrite(baseline,
			Frag(300f, 0.5f, ProductType.b, 3, 1, ProductType.y, 5, -18.0, 2, excludeFromQuant: true)), Is.Negative);
	}
}
