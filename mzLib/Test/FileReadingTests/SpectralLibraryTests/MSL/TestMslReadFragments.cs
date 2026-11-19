using MassSpectrometry;
using Microsoft.Win32.SafeHandles;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;
using System.IO.MemoryMappedFiles;
using System.Linq;
using System.Reflection;
using System.Threading.Tasks;

namespace Test.MslSpectralLibrary;

/// <summary>
/// Tests for <see cref="MslLibrary.GetFragmentCount"/> and <see cref="MslLibrary.ReadFragments"/>:
/// the same fragments as <see cref="MslLibrary.GetEntry(int)"/>, copied into a caller-owned buffer
/// without allocating, in every load mode.
/// </summary>
[TestFixture]
public sealed class TestMslReadFragments
{
	public enum Mode { FullLoad, IndexOnly, CompressedFallback }

	private const int OrdinaryEntries = 20;
	private const int ManyFragments = 150; // more than two of the reader's 64-record chunks
	private const double CustomLossA = 123.4567;
	private const double CustomLossB = 98.7654;

	private string _dir = null!;
	private string _path = null!;
	private string _compressedPath = null!;

	[OneTimeSetUp]
	public void OneTimeSetUp()
	{
		_dir = Path.Combine(Path.GetTempPath(), "MslReadFragments_" + Guid.NewGuid().ToString("N"));
		Directory.CreateDirectory(_dir);
		_path = Path.Combine(_dir, "lib.msl");
		_compressedPath = Path.Combine(_dir, "lib_zstd.msl");
		MslWriter.Write(_path, BuildEntries());
		MslWriter.Write(_compressedPath, BuildEntries(), compressionLevel: 3);
	}

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		try { Directory.Delete(_dir, recursive: true); } catch (IOException) { }
	}

	private static MslFragmentIon Frag(float mz, ProductType type, int number, int charge = 1,
		ProductType? secondaryType = null, int secondaryNumber = 0, double neutralLoss = 0.0,
		int residuePosition = 0, bool excludeFromQuant = false, float intensity = 0.5f) => new()
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

	private static MslLibraryEntry Entry(int i, List<MslFragmentIon> fragments)
	{
		string seq = "PEPTIDE" + new string('K', i % 5 + 1) + (char)('A' + i % 20);
		return new MslLibraryEntry
		{
			FullSequence = seq,
			BaseSequence = seq,
			PrecursorMz = 400.0 + i * 7.25,
			ChargeState = 2 + i % 2,
			RetentionTime = 10.0 + i,
			ProteinAccession = "P" + i.ToString("D5"),
			DissociationType = DissociationType.HCD,
			Nce = 28,
			MatchedFragmentIons = fragments
		};
	}

	/// <summary>
	/// Entry 0 has every kind of fragment, including two custom losses; entry 1 has none; entry 2
	/// has more fragments than one read chunk; the rest are ordinary.
	/// </summary>
	private static List<MslLibraryEntry> BuildEntries()
	{
		var entries = new List<MslLibraryEntry>
		{
			Entry(0, new List<MslFragmentIon>
			{
				Frag(138.06f, ProductType.D, 0),                                              // diagnostic
				Frag(200.10f, ProductType.b, 3, residuePosition: 3),
				Frag(250.10f, ProductType.y, 4, neutralLoss: 18.010565, residuePosition: 5),  // H2O
				Frag(251.10f, ProductType.y, 4, neutralLoss: 17.026549, residuePosition: 5),  // NH3
				Frag(300.20f, ProductType.y, 5, residuePosition: 4, intensity: 1.0f),
				Frag(350.10f, ProductType.b, 5, neutralLoss: 97.976895, residuePosition: 5),  // H3PO4
				Frag(370.10f, ProductType.b, 6, neutralLoss: CustomLossA, residuePosition: 6),
				Frag(371.10f, ProductType.y, 6, neutralLoss: CustomLossB, residuePosition: 3),
				Frag(380.10f, ProductType.b, 7, excludeFromQuant: true),
				Frag(390.10f, ProductType.y, 7, charge: 2),
				Frag(400.30f, ProductType.b, 3, secondaryType: ProductType.y, secondaryNumber: 6), // internal
			}),
			Entry(1, new List<MslFragmentIon>()),
			Entry(2, Enumerable.Range(0, ManyFragments)
				.Select(k => Frag(150f + 3.5f * k, k % 2 == 0 ? ProductType.b : ProductType.y, k / 2 + 1,
					residuePosition: k % 9, intensity: 0.1f + k % 7 * 0.1f))
				.ToList()),
		};

		for (int i = 3; i < 3 + OrdinaryEntries; i++)
			entries.Add(Entry(i, new List<MslFragmentIon>
			{
				Frag(200f + i, ProductType.b, 2, residuePosition: 2),
				Frag(300f + i, ProductType.y, 3, residuePosition: 6, intensity: 1.0f),
				Frag(400f + i, ProductType.y, 4, neutralLoss: 18.010565, residuePosition: 5),
			}));

		return entries;
	}

	private MslLibrary Open(Mode mode) => mode switch
	{
		Mode.FullLoad => MslLibrary.Load(_path),
		Mode.IndexOnly => MslLibrary.LoadIndexOnly(_path),
		Mode.CompressedFallback => MslLibrary.LoadIndexOnly(_compressedPath),
		_ => throw new ArgumentOutOfRangeException(nameof(mode))
	};

	/// <summary>
	/// The stored fields of <see cref="MslFragmentIon"/>: its settable properties. The get-only ones
	/// (<c>IsInternalFragment</c>, <c>IsDiagnosticIon</c>, <c>Annotation</c>) are computed from them.
	/// </summary>
	private static readonly PropertyInfo[] StoredFields = typeof(MslFragmentIon)
		.GetProperties(BindingFlags.Public | BindingFlags.Instance)
		.Where(p => p.CanWrite)
		.ToArray();

	/// <summary>Asserts every stored field of <paramref name="ion"/> equals the value's.</summary>
	private static void AssertSameFragment(MslFragmentValue value, MslFragmentIon ion, string where)
	{
		Assert.That(value.IsInternalFragment, Is.EqualTo(ion.IsInternalFragment), $"{where}: IsInternalFragment");
		Assert.That(value.IsDiagnosticIon, Is.EqualTo(ion.IsDiagnosticIon), $"{where}: IsDiagnosticIon");
		foreach (PropertyInfo p in StoredFields)
		{
			object? expected = p.GetValue(ion);
			object? actual = typeof(MslFragmentValue).GetProperty(p.Name)!.GetValue(value);
			Assert.That(actual, Is.EqualTo(expected), $"{where}: {p.Name}");
		}
	}

	// ── Same fragments as GetEntry ─────────────────────────────────────────────

	[TestCase(Mode.FullLoad)]
	[TestCase(Mode.IndexOnly)]
	[TestCase(Mode.CompressedFallback)]
	public void ReadFragments_GivesTheFragmentsOfGetEntry_ForEveryEntry(Mode mode)
	{
		using MslLibrary lib = Open(mode);
		Assert.That(lib.IsIndexOnly, Is.EqualTo(mode == Mode.IndexOnly));

		var buffer = new MslFragmentValue[ManyFragments + 5];
		for (int i = 0; i < lib.PrecursorCount; i++)
		{
			List<MslFragmentIon> ions = lib.GetEntry(i)!.MatchedFragmentIons;

			Assert.That(lib.GetFragmentCount(i), Is.EqualTo(ions.Count), $"count of entry {i}");
			Assert.That(lib.ReadFragments(i, buffer), Is.EqualTo(ions.Count), $"written for entry {i}");
			for (int k = 0; k < ions.Count; k++)
				AssertSameFragment(buffer[k], ions[k], $"entry {i} fragment {k}");
		}
	}

	[TestCase(Mode.FullLoad)]
	[TestCase(Mode.IndexOnly)]
	[TestCase(Mode.CompressedFallback)]
	public void ReadFragments_DecodesNamedAndCustomLosses(Mode mode)
	{
		using MslLibrary lib = Open(mode);
		var buffer = new MslFragmentValue[lib.GetFragmentCount(0)];
		lib.ReadFragments(0, buffer);

		double[] losses = buffer.Select(f => f.NeutralLoss).Where(l => l != 0).OrderBy(l => l).ToArray();
		Assert.That(losses, Is.EqualTo(new[] { 17.026549, 18.010565, 97.976895, CustomLossB, CustomLossA }).Within(1e-4));

		// A custom loss uses the residue field on disk for its table index, so it reads back as 0
		Assert.That(buffer.Where(f => Math.Abs(f.NeutralLoss - CustomLossA) < 1e-4 || Math.Abs(f.NeutralLoss - CustomLossB) < 1e-4)
			.Select(f => f.ResiduePosition), Is.All.Zero);
		Assert.That(buffer.Count(f => f.IsInternalFragment), Is.EqualTo(1));
		Assert.That(buffer.Count(f => f.IsDiagnosticIon), Is.EqualTo(1));
		Assert.That(buffer.Count(f => f.ExcludeFromQuant), Is.EqualTo(1));
		Assert.That(buffer.Count(f => f.Charge == 2), Is.EqualTo(1));
	}

	/// <summary>
	/// The library-level test above reaches the memory map on a local drive. Both index-only read
	/// paths are checked here against the list-returning reader.
	/// </summary>
	[TestCase(true)]
	[TestCase(false)]
	public void IndexOnlyReadPaths_MatchLoadFragmentsOnDemand(bool memoryMapFragments)
	{
		using MslLibraryData raw = MslReader.LoadIndexOnly(_path, memoryMapFragments);
		var buffer = new MslFragmentValue[ManyFragments];

		for (int i = 0; i < raw.Count; i++)
		{
			List<MslFragmentIon> ions = raw.LoadFragmentsOnDemand(i);
			Assert.That(raw.GetFragmentCount(i), Is.EqualTo(ions.Count));
			Assert.That(raw.ReadFragments(i, buffer), Is.EqualTo(ions.Count));
			for (int k = 0; k < ions.Count; k++)
				AssertSameFragment(buffer[k], ions[k], $"entry {i} fragment {k}");
		}
	}

	// ── No allocation, no cache ────────────────────────────────────────────────

	[TestCase(Mode.FullLoad)]
	[TestCase(Mode.IndexOnly)]
	[TestCase(Mode.CompressedFallback)]
	public void ReadFragments_AllocatesNothing(Mode mode)
	{
		using MslLibrary lib = Open(mode);
		var buffer = new MslFragmentValue[ManyFragments];

		long fragments = ReadAll(lib, buffer); // warm up
		long before = GC.GetAllocatedBytesForCurrentThread();
		for (int pass = 0; pass < 10; pass++)
			fragments += ReadAll(lib, buffer);
		long allocated = GC.GetAllocatedBytesForCurrentThread() - before;

		Assert.That(fragments, Is.GreaterThan(0));
		Assert.That(allocated, Is.Zero);

		static long ReadAll(MslLibrary lib, Span<MslFragmentValue> buffer)
		{
			long n = 0;
			for (int i = 0; i < lib.PrecursorCount; i++)
				n += lib.GetFragmentCount(i) + lib.ReadFragments(i, buffer);
			return n;
		}
	}

	[TestCase(true)]
	[TestCase(false)]
	public void IndexOnlyReadPaths_AllocateNothing(bool memoryMapFragments)
	{
		using MslLibraryData raw = MslReader.LoadIndexOnly(_path, memoryMapFragments);
		var buffer = new MslFragmentValue[ManyFragments];

		for (int i = 0; i < raw.Count; i++)
			raw.ReadFragments(i, buffer); // warm up

		long before = GC.GetAllocatedBytesForCurrentThread();
		for (int pass = 0; pass < 10; pass++)
			for (int i = 0; i < raw.Count; i++)
				raw.ReadFragments(i, buffer);

		Assert.That(GC.GetAllocatedBytesForCurrentThread() - before, Is.Zero);
	}

	[Test]
	public void IndexOnly_ReadFragments_DoesNotTouchTheEntryCache()
	{
		using MslLibrary lib = MslLibrary.LoadIndexOnly(_path);
		var buffer = new MslFragmentValue[ManyFragments];

		for (int i = 0; i < lib.PrecursorCount; i++)
		{
			lib.GetFragmentCount(i);
			lib.ReadFragments(i, buffer);
		}

		Assert.That(lib.IndexStatisticsForTesting.LruHits, Is.Zero);
		Assert.That(lib.IndexStatisticsForTesting.LruMisses, Is.Zero);
	}

	[Test]
	public void IndexOnly_ParallelReads_MatchSequentialReads()
	{
		using MslLibrary lib = MslLibrary.LoadIndexOnly(_path);
		int n = lib.PrecursorCount;

		var expected = new MslFragmentValue[n][];
		for (int i = 0; i < n; i++)
		{
			expected[i] = new MslFragmentValue[lib.GetFragmentCount(i)];
			lib.ReadFragments(i, expected[i]);
		}

		var actual = new MslFragmentValue[n * 40][];
		Parallel.For(0, actual.Length, new ParallelOptions { MaxDegreeOfParallelism = 8 }, j =>
		{
			int i = j % n;
			var buffer = new MslFragmentValue[lib.GetFragmentCount(i)];
			lib.ReadFragments(i, buffer);
			actual[j] = buffer;
		});

		for (int j = 0; j < actual.Length; j++)
			Assert.That(actual[j], Is.EqualTo(expected[j % n]), $"read {j}");
	}

	// ── Failure contract ───────────────────────────────────────────────────────

	[TestCase(Mode.FullLoad)]
	[TestCase(Mode.IndexOnly)]
	[TestCase(Mode.CompressedFallback)]
	public void EntryWithoutFragments_ReadsZero(Mode mode)
	{
		using MslLibrary lib = Open(mode);
		Assert.That(lib.GetFragmentCount(1), Is.Zero);
		Assert.That(lib.ReadFragments(1, Span<MslFragmentValue>.Empty), Is.Zero);
	}

	[TestCase(Mode.FullLoad)]
	[TestCase(Mode.IndexOnly)]
	[TestCase(Mode.CompressedFallback)]
	public void ShortDestination_Throws_AndWritesNothing(Mode mode)
	{
		using MslLibrary lib = Open(mode);
		int count = lib.GetFragmentCount(0);
		var buffer = new MslFragmentValue[count - 1];

		Assert.That(() => lib.ReadFragments(0, buffer),
			Throws.ArgumentException.With.Property("ParamName").EqualTo("destination"));
		Assert.That(buffer, Is.All.EqualTo(default(MslFragmentValue)));
	}

	[TestCase(Mode.FullLoad)]
	[TestCase(Mode.IndexOnly)]
	[TestCase(Mode.CompressedFallback)]
	public void OutOfRangeIndex_Throws(Mode mode)
	{
		using MslLibrary lib = Open(mode);
		var buffer = new MslFragmentValue[ManyFragments];

		foreach (int bad in new[] { -1, lib.PrecursorCount })
		{
			Assert.That(() => lib.GetFragmentCount(bad),
				Throws.TypeOf<ArgumentOutOfRangeException>().With.Property("ParamName").EqualTo("precursorIdx"));
			Assert.That(() => lib.ReadFragments(bad, buffer),
				Throws.TypeOf<ArgumentOutOfRangeException>().With.Property("ParamName").EqualTo("precursorIdx"));
		}
	}

	[TestCase(Mode.FullLoad)]
	[TestCase(Mode.IndexOnly)]
	[TestCase(Mode.CompressedFallback)]
	public void AfterDispose_Throws(Mode mode)
	{
		MslLibrary lib = Open(mode);
		lib.Dispose();
		var buffer = new MslFragmentValue[ManyFragments];

		Assert.That(() => lib.GetFragmentCount(0), Throws.TypeOf<ObjectDisposedException>());
		Assert.That(() => lib.ReadFragments(0, buffer), Throws.TypeOf<ObjectDisposedException>());
	}

	[Test]
	public void FullLoadData_RefusesFileReads()
	{
		MslLibraryData raw = MslReader.Load(_path);
		var buffer = new MslFragmentValue[ManyFragments];

		Assert.That(() => raw.GetFragmentCount(0), Throws.InvalidOperationException);
		Assert.That(() => raw.ReadFragments(0, buffer), Throws.InvalidOperationException);
	}

	// ── MslFragmentValue ───────────────────────────────────────────────────────

	/// <summary>
	/// The value carries every stored field of <see cref="MslFragmentIon"/>, with the same type, so
	/// a field added to the class cannot be left out of the value without this test failing.
	/// </summary>
	[Test]
	public void FragmentValue_HasEveryStoredFieldOfFragmentIon()
	{
		Dictionary<string, Type> valueProperties = typeof(MslFragmentValue)
			.GetProperties(BindingFlags.Public | BindingFlags.Instance)
			.ToDictionary(p => p.Name, p => p.PropertyType);

		Assert.That(StoredFields, Has.Length.GreaterThanOrEqualTo(10)); // the ten fields at the time of writing
		foreach (PropertyInfo p in StoredFields)
		{
			Assert.That(valueProperties.TryGetValue(p.Name, out Type? type), Is.True, $"missing {p.Name}");
			Assert.That(type, Is.EqualTo(p.PropertyType), p.Name);
		}
	}

	[Test]
	public void FragmentValue_FromIon_CopiesEveryField()
	{
		MslFragmentIon ion = Frag(512.25f, ProductType.c, 7, charge: 3, secondaryType: ProductType.zDot,
			secondaryNumber: 11, neutralLoss: 44.5, residuePosition: 9, excludeFromQuant: true, intensity: 0.75f);

		AssertSameFragment(new MslFragmentValue(ion), ion, "from ion");
		Assert.That(() => new MslFragmentValue(null!), Throws.ArgumentNullException);
	}

	// ── Guards below MslLibrary ────────────────────────────────────────────────
	// MslLibrary refuses a disposed or out-of-range read before the data layer runs, so these
	// tests call the data layer and the record readers directly.

	[TestCase(true)]
	[TestCase(false)]
	public void DataLayer_RefusesReadsAfterDispose(bool memoryMapFragments)
	{
		MslLibraryData raw = MslReader.LoadIndexOnly(_path, memoryMapFragments);
		raw.Dispose();
		var buffer = new MslFragmentValue[ManyFragments];

		Assert.That(() => raw.ReadFragments(0, buffer), Throws.TypeOf<ObjectDisposedException>());
	}

	[TestCase(true)]
	[TestCase(false)]
	public void DataLayer_RefusesAnOutOfRangeIndex(bool memoryMapFragments)
	{
		using MslLibraryData raw = MslReader.LoadIndexOnly(_path, memoryMapFragments);
		var buffer = new MslFragmentValue[ManyFragments];

		foreach (int bad in new[] { -1, raw.Count })
		{
			Assert.That(() => raw.GetFragmentCount(bad), Throws.TypeOf<ArgumentOutOfRangeException>());
			Assert.That(() => raw.ReadFragments(bad, buffer), Throws.TypeOf<ArgumentOutOfRangeException>());
		}
	}

	[Test]
	public void MappedRecordRead_PastTheEndOfTheFile_Throws()
	{
		long length = new FileInfo(_path).Length;
		using var map = MemoryMappedFile.CreateFromFile(_path, FileMode.Open, null, 0, MemoryMappedFileAccess.Read);
		using MemoryMappedViewAccessor view = map.CreateViewAccessor(0, 0, MemoryMappedFileAccess.Read);
		var records = new MslFragmentRecord[2];

		Assert.That(() => MslReader.ReadFragmentRecordsAt(view, length, length - MslFormat.FragmentRecordSize, records),
			Throws.TypeOf<EndOfStreamException>());
		Assert.That(() => MslReader.ReadFragmentRecordsAt(view, length, -1, records),
			Throws.TypeOf<EndOfStreamException>());
	}

	[Test]
	public void PositionalRecordRead_PastTheEndOfTheFile_Throws()
	{
		long length = new FileInfo(_path).Length;
		using SafeFileHandle handle = File.OpenHandle(_path, FileMode.Open, FileAccess.Read, FileShare.Read);
		var records = new MslFragmentRecord[2];

		// One record is left before the end: the first read returns it, the next returns nothing
		Assert.That(() => MslReader.ReadFragmentRecordsAt(handle, length - MslFormat.FragmentRecordSize, records),
			Throws.TypeOf<EndOfStreamException>());
	}

	/// <summary>
	/// A custom loss whose table index is missing from the file's table decodes as no loss, as
	/// GetEntry's decoder does (a defensive fallback; a well-formed file never has one).
	/// </summary>
	[TestCase((short)0)]
	[TestCase((short)99)]
	public void CustomLossIndexOutsideTheTable_DecodesAsNoLoss(short tableIndex)
	{
		var record = new MslFragmentRecord
		{
			Mz = 300f,
			Intensity = 1f,
			ProductType = (short)ProductType.y,
			SecondaryProductType = -1,
			FragmentNumber = 3,
			ResiduePosition = tableIndex,
			Charge = 1,
			Flags = MslFormat.EncodeFragmentFlags(false, false, MslFormat.NeutralLossCode.Custom, false)
		};

		MslFragmentValue value = MslReader.ToFragmentValue(in record, new[] { 0.0, CustomLossA });

		Assert.That(value.NeutralLoss, Is.Zero);
		Assert.That(value.ResiduePosition, Is.Zero);
	}

	[Test]
	public void FullLoadEntryWithNoFragmentList_ReadsZero()
	{
		using MslLibrary lib = MslLibrary.Load(_path);
		lib.GetEntry(3)!.MatchedFragmentIons = null!; // full load returns the stored entry itself

		Assert.That(lib.GetFragmentCount(3), Is.Zero);
		Assert.That(lib.ReadFragments(3, Span<MslFragmentValue>.Empty), Is.Zero);
	}
}
