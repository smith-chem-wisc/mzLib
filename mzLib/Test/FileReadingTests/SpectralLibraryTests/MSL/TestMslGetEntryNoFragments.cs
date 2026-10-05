using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;

namespace Test.MslSpectralLibrary;

/// <summary>
/// Tests for <see cref="MslLibrary.GetEntry(int, bool)"/> and for the index-only loader no
/// longer storing fragments on the library's shared entry skeletons.
/// </summary>
[TestFixture]
public sealed class TestMslGetEntryNoFragments
{
	private const int EntryCount = 40;

	private string _dir = null!;
	private string _path = null!;

	[OneTimeSetUp]
	public void OneTimeSetUp()
	{
		_dir = Path.Combine(Path.GetTempPath(), "MslGetEntryNoFragments_" + Guid.NewGuid().ToString("N"));
		Directory.CreateDirectory(_dir);
		_path = Path.Combine(_dir, "lib.msl");
		MslLibrary.Save(_path, BuildEntries());
	}

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		try { Directory.Delete(_dir, recursive: true); } catch (IOException) { }
	}

	private static List<MslLibraryEntry> BuildEntries()
	{
		var entries = new List<MslLibraryEntry>();
		for (int i = 0; i < EntryCount; i++)
		{
			string seq = "PEPTIDE" + new string('K', i % 5 + 1) + (char)('A' + i % 20);
			entries.Add(new MslLibraryEntry
			{
				FullSequence = seq,
				BaseSequence = seq,
				PrecursorMz = 400.0 + i * 7.25,
				ChargeState = 2 + i % 2,
				RetentionTime = 10.0 + i,
				IsDecoy = i % 4 == 3,
				ProteinAccession = "P" + i.ToString("D5"),
				ProteinName = "Protein " + i,
				GeneName = "GENE" + i,
				DissociationType = DissociationType.HCD,
				Nce = 28,
				MatchedFragmentIons = new List<MslFragmentIon>
				{
					new MslFragmentIon { Mz = 100f + i, Intensity = 1000f, ProductType = ProductType.b, FragmentNumber = 2, Charge = 1 },
					new MslFragmentIon { Mz = 200f + i, Intensity = 500f, ProductType = ProductType.y, FragmentNumber = 3, Charge = 1 }
				}
			});
		}
		return entries;
	}

	private static long LruMisses(MslLibrary lib)
	{
		var index = (MslIndex)typeof(MslLibrary)
			.GetField("_index", BindingFlags.NonPublic | BindingFlags.Instance)!
			.GetValue(lib)!;
		return index.GetStatistics().LruMisses;
	}

	private static void AssertSameMetadata(MslLibraryEntry expected, MslLibraryEntry actual)
	{
		Assert.That(actual.FullSequence, Is.EqualTo(expected.FullSequence));
		Assert.That(actual.ChargeState, Is.EqualTo(expected.ChargeState));
		Assert.That(actual.PrecursorMz, Is.EqualTo(expected.PrecursorMz));
		Assert.That(actual.RetentionTime, Is.EqualTo(expected.RetentionTime));
		Assert.That(actual.IsDecoy, Is.EqualTo(expected.IsDecoy));
		Assert.That(actual.ProteinAccession, Is.EqualTo(expected.ProteinAccession));
		Assert.That(actual.ProteinName, Is.EqualTo(expected.ProteinName));
		Assert.That(actual.GeneName, Is.EqualTo(expected.GeneName));
	}

	[Test]
	public void IndexOnly_NoFragments_ReturnsMetadataWithoutReadingFragments()
	{
		using MslLibrary full = MslLibrary.Load(_path);
		using MslLibrary indexOnly = MslLibrary.LoadIndexOnly(_path);
		Assert.That(indexOnly.IsIndexOnly, Is.True);

		long missesBefore = LruMisses(indexOnly);
		for (int i = 0; i < EntryCount; i++)
		{
			MslLibraryEntry? meta = indexOnly.GetEntry(i, includeFragments: false);
			Assert.That(meta, Is.Not.Null);
			AssertSameMetadata(full.GetEntry(i)!, meta!);
			Assert.That(meta!.MatchedFragmentIons, Is.Empty);
		}

		Assert.That(LruMisses(indexOnly), Is.EqualTo(missesBefore),
			"includeFragments: false must not go through the entry loader.");
	}

	[Test]
	public void IndexOnly_NoFragments_StaysEmptyAfterFullReadOfSameEntry()
	{
		using MslLibrary indexOnly = MslLibrary.LoadIndexOnly(_path);

		for (int i = 0; i < EntryCount; i++)
		{
			MslLibraryEntry withFragments = indexOnly.GetEntry(i)!;
			Assert.That(withFragments.MatchedFragmentIons, Has.Count.EqualTo(2));

			MslLibraryEntry meta = indexOnly.GetEntry(i, includeFragments: false)!;
			Assert.That(meta, Is.Not.SameAs(withFragments));
			Assert.That(meta.MatchedFragmentIons, Is.Empty,
				"A fragment read must not leave its fragments on the shared skeleton.");
			AssertSameMetadata(withFragments, meta);
		}
	}

	[Test]
	public void IndexOnly_FragmentReads_MatchFullLoad()
	{
		using MslLibrary full = MslLibrary.Load(_path);
		using MslLibrary indexOnly = MslLibrary.LoadIndexOnly(_path);

		for (int pass = 0; pass < 2; pass++)
			for (int i = 0; i < EntryCount; i++)
			{
				List<MslFragmentIon> expected = full.GetEntry(i)!.MatchedFragmentIons;
				List<MslFragmentIon> actual = indexOnly.GetEntry(i, includeFragments: true)!.MatchedFragmentIons;
				Assert.That(actual.Select(f => f.Mz), Is.EqualTo(expected.Select(f => f.Mz)));
				Assert.That(actual.Select(f => f.Intensity), Is.EqualTo(expected.Select(f => f.Intensity)));
			}
	}

	[Test]
	public void FullLoad_NoFragments_ReturnsSameEntryAsGetEntry()
	{
		using MslLibrary full = MslLibrary.Load(_path);

		for (int i = 0; i < EntryCount; i++)
			Assert.That(full.GetEntry(i, includeFragments: false), Is.SameAs(full.GetEntry(i)));
	}

	[TestCase(false)]
	[TestCase(true)]
	public void NoFragments_OutOfRange_ReturnsNull(bool indexOnlyMode)
	{
		using MslLibrary lib = indexOnlyMode ? MslLibrary.LoadIndexOnly(_path) : MslLibrary.Load(_path);

		Assert.That(lib.GetEntry(-1, includeFragments: false), Is.Null);
		Assert.That(lib.GetEntry(EntryCount, includeFragments: false), Is.Null);
	}

	[Test]
	public void NoFragments_AfterDispose_Throws()
	{
		MslLibrary lib = MslLibrary.LoadIndexOnly(_path);
		lib.Dispose();

		Assert.Throws<ObjectDisposedException>(() => lib.GetEntry(0, includeFragments: false));
	}

	[Test]
	public void CalibratedIndexOnly_NoFragments_ReturnsMetadata()
	{
		using MslLibrary indexOnly = MslLibrary.LoadIndexOnly(_path);
		using MslLibrary calibrated = indexOnly.WithCalibratedRetentionTimes(2.0, 1.0);

		MslLibraryEntry meta = calibrated.GetEntry(5, includeFragments: false)!;
		Assert.That(meta.ProteinAccession, Is.EqualTo(indexOnly.GetEntry(5)!.ProteinAccession));
		Assert.That(meta.MatchedFragmentIons, Is.Empty);
	}

	/// <summary>
	/// Guards the property list in <c>CopyWithFragments</c>: a property added to
	/// <see cref="MslLibraryEntry"/> later must be copied too, or index-only entries lose it.
	/// </summary>
	[Test]
	public void CopyWithFragments_CopiesEveryWritableProperty()
	{
		var source = new MslLibraryEntry();
		PropertyInfo[] props = typeof(MslLibraryEntry)
			.GetProperties(BindingFlags.Public | BindingFlags.Instance)
			.Where(p => p.CanRead && p.CanWrite && p.Name != nameof(MslLibraryEntry.MatchedFragmentIons))
			.ToArray();

		foreach (PropertyInfo p in props)
			p.SetValue(source, NonDefaultValue(p));

		var fragments = new List<MslFragmentIon> { new MslFragmentIon { Mz = 1f } };
		MslLibraryEntry copy = MslLibrary.CopyWithFragments(source, fragments);

		Assert.That(copy.MatchedFragmentIons, Is.SameAs(fragments));
		foreach (PropertyInfo p in props)
			Assert.That(p.GetValue(copy), Is.EqualTo(p.GetValue(source)), $"{p.Name} was not copied.");
	}

	private static object NonDefaultValue(PropertyInfo p)
	{
		Type t = p.PropertyType;
		if (t == typeof(string)) return "value of " + p.Name;
		if (t == typeof(int)) return 7;
		if (t == typeof(double)) return 3.5;
		if (t == typeof(float)) return 2.5f;
		if (t == typeof(bool)) return true;
		if (t.IsEnum)
		{
			Array values = Enum.GetValues(t);
			return values.GetValue(values.Length - 1)!;
		}
		throw new NotSupportedException(
			$"Add a non-default value for {t.Name} ({p.Name}) to this test.");
	}
}
