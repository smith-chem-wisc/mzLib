using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using ZstdSharp;

namespace Test.MslSpectralLibrary;

/// <summary>
/// Tests for reading and writing .msl files in bounded pieces rather than as one buffer:
/// <see cref="MslReader.Load"/> streams the file, the writer compresses through a stream, and
/// index-only fragment reads are positional (no shared file position, no lock).
/// None of these paths needs a file over 2 GB to be exercised; the 2 GB limit they remove came
/// from holding the whole file (or the whole fragment section) in one <c>byte[]</c>.
/// </summary>
[TestFixture]
public sealed class TestMslStreamedReads
{
	private static readonly string OutputDirectory =
		Path.Combine(Path.GetTempPath(), "MslStreamedReadsTests");

	[OneTimeSetUp]
	public void OneTimeSetUp() => Directory.CreateDirectory(OutputDirectory);

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(OutputDirectory))
			Directory.Delete(OutputDirectory, recursive: true);
	}

	private static string TempPath(string name) => Path.Combine(OutputDirectory, name + ".msl");

	/// <summary>
	/// Entries whose fragment m/z values are unique to the entry, so a fragment list attached
	/// to the wrong precursor is detectable.
	/// </summary>
	private static List<MslLibraryEntry> MakeEntries(int count)
	{
		var entries = new List<MslLibraryEntry>(count);
		for (int i = 0; i < count; i++)
		{
			string seq = "PEPTIDE" + new string('K', i % 7) + "R" + i;
			int nFragments = 1 + i % 9;
			entries.Add(new MslLibraryEntry
			{
				FullSequence = seq,
				BaseSequence = seq,
				PrecursorMz = 400.0 + i * 0.37,
				ChargeState = 2 + i % 2,
				RetentionTime = 10.0 + i * 0.05,
				MoleculeType = MslFormat.MoleculeType.Peptide,
				DissociationType = DissociationType.HCD,
				MatchedFragmentIons = Enumerable.Range(0, nFragments).Select(f => new MslFragmentIon
				{
					Mz = 100f + i + f * 0.001f,
					Intensity = 1f - f * 0.01f,
					ProductType = ProductType.y,
					FragmentNumber = f + 1,
					ResiduePosition = f + 1,
					Charge = 1
				}).ToList()
			});
		}
		return entries;
	}

	private static float[] FragmentMzs(MslLibraryEntry entry) =>
		entry.MatchedFragmentIons.Select(f => (float)f.Mz).ToArray();

	private static void AssertSameFragments(IReadOnlyList<MslLibraryEntry> expected, IReadOnlyList<MslLibraryEntry> actual)
	{
		Assert.That(actual.Count, Is.EqualTo(expected.Count));
		var bySequence = expected.ToDictionary(e => e.FullSequence);
		foreach (MslLibraryEntry entry in actual)
			Assert.That(FragmentMzs(entry), Is.EqualTo(FragmentMzs(bySequence[entry.FullSequence])),
				$"fragments of {entry.FullSequence}");
	}

	/// <summary>
	/// The writer lays fragment blocks out in precursor order, so a full load is one forward
	/// pass. A file whose blocks are not in that order must still attach each block to its own
	/// precursor: swap two precursor records (so their offsets descend) and re-seal the CRC.
	/// </summary>
	[TestCase(0)]
	[TestCase(3)]
	public void Load_FragmentBlocksOutOfPrecursorOrder_EachEntryKeepsItsOwnFragments(int compressionLevel)
	{
		string path = TempPath(nameof(Load_FragmentBlocksOutOfPrecursorOrder_EachEntryKeepsItsOwnFragments) + compressionLevel);
		List<MslLibraryEntry> written = MakeEntries(5);
		MslWriter.Write(path, written, compressionLevel);

		long precursorSection = MslReader.ReadHeaderOnly(path).PrecursorSectionOffset;
		byte[] bytes = File.ReadAllBytes(path);
		int first = (int)precursorSection;
		int last = first + 4 * MslFormat.PrecursorRecordSize;
		byte[] firstRecord = bytes.AsSpan(first, MslFormat.PrecursorRecordSize).ToArray();
		bytes.AsSpan(last, MslFormat.PrecursorRecordSize).CopyTo(bytes.AsSpan(first));
		firstRecord.CopyTo(bytes.AsSpan(last));
		RecomputeCrc(bytes);
		File.WriteAllBytes(path, bytes);

		MslLibraryData loaded = MslReader.Load(path);

		Assert.That(loaded.Entries[0].FullSequence, Is.EqualTo(written[4].FullSequence), "records were swapped");
		AssertSameFragments(written, loaded.Entries);
	}

	/// <summary>
	/// Both writers now compress through a stream. The frame must still decode with the
	/// whole-buffer <see cref="Decompressor.Unwrap(ReadOnlySpan{byte}, Span{byte})"/> call that
	/// readers before this change use, and must decode to exactly the uncompressed size.
	/// </summary>
	[TestCase(false)]
	[TestCase(true)]
	public void CompressedWrite_FrameDecodesWithWholeBufferDecoder(bool streamingWriter)
	{
		string path = TempPath(nameof(CompressedWrite_FrameDecodesWithWholeBufferDecoder) + streamingWriter);
		List<MslLibraryEntry> written = MakeEntries(200);
		if (streamingWriter)
			MslWriter.WriteStreaming(path, written, compressionLevel: 3);
		else
			MslWriter.Write(path, written, compressionLevel: 3);

		byte[] bytes = File.ReadAllBytes(path);
		long compressedSize = BitConverter.ToInt64(bytes, MslFormat.HeaderSize);
		long uncompressedSize = BitConverter.ToInt64(bytes, MslFormat.HeaderSize + 8);
		long frameStart = MslReader.ReadHeaderOnly(path).FragmentSectionOffset;

		Assert.That(uncompressedSize, Is.EqualTo(written.Sum(e => e.MatchedFragmentIons.Count) * MslFormat.FragmentRecordSize));

		using var decompressor = new Decompressor();
		byte[] decompressed = new byte[uncompressedSize];
		int decoded = decompressor.Unwrap(bytes.AsSpan((int)frameStart, (int)compressedSize), decompressed);
		Assert.That(decoded, Is.EqualTo(uncompressedSize));

		AssertSameFragments(written, MslReader.Load(path).Entries);
	}

	/// <summary>
	/// Index-only fragment reads are positional, so many threads can read at once. Every
	/// entry fetched concurrently must match the full load.
	/// </summary>
	[Test]
	public void LoadIndexOnly_ParallelGetEntry_MatchesFullLoad()
	{
		string path = TempPath(nameof(LoadIndexOnly_ParallelGetEntry_MatchesFullLoad));
		List<MslLibraryEntry> written = MakeEntries(500);
		MslWriter.Write(path, written);

		using MslLibrary full = MslLibrary.Load(path);
		using MslLibrary indexOnly = MslLibrary.LoadIndexOnly(path);
		Assert.That(indexOnly.IsIndexOnly, Is.True);

		int mismatches = 0;
		Parallel.For(0, 4 * written.Count, new ParallelOptions { MaxDegreeOfParallelism = 8 }, k =>
		{
			int i = k % written.Count;
			if (!FragmentMzs(indexOnly.GetEntry(i)!).SequenceEqual(FragmentMzs(full.GetEntry(i)!)))
				Interlocked.Increment(ref mismatches);
		});

		Assert.That(mismatches, Is.EqualTo(0));
	}

	/// <summary>
	/// A compressed file whose frame is cut short must fail to load rather than return
	/// partial fragments. The CRC is re-sealed so the truncation reaches the decoder.
	/// </summary>
	[Test]
	public void Load_CompressedFrameShorterThanDescriptorSays_Throws()
	{
		string path = TempPath(nameof(Load_CompressedFrameShorterThanDescriptorSays_Throws));
		MslWriter.Write(path, MakeEntries(50), compressionLevel: 3);

		byte[] bytes = File.ReadAllBytes(path);
		long compressedSize = BitConverter.ToInt64(bytes, MslFormat.HeaderSize);
		BitConverter.GetBytes(compressedSize / 2).CopyTo(bytes, MslFormat.HeaderSize);
		RecomputeCrc(bytes);
		File.WriteAllBytes(path, bytes);

		Assert.That(() => MslReader.Load(path), Throws.InstanceOf<Exception>());
	}

	/// <summary>
	/// The slicing-by-8 CRC must equal CRC-32/ISO-HDLC exactly: the standard check value, and
	/// a bit-at-a-time reference over every length 0–64 plus a large buffer, fed whole or in
	/// uneven pieces. Files written before this change carry CRCs from the byte-table version.
	/// </summary>
	[Test]
	public void MslCrc32_MatchesBitwiseReference()
	{
		Assert.That(MslCrc32.Compute("123456789"u8), Is.EqualTo(0xCBF43926u));

		var random = new Random(42);
		byte[] data = new byte[1 << 20];
		random.NextBytes(data);

		for (int length = 0; length <= 64; length++)
			Assert.That(MslCrc32.Compute(data.AsSpan(0, length)), Is.EqualTo(BitwiseCrc32(data.AsSpan(0, length))), $"length {length}");

		uint expected = BitwiseCrc32(data);
		Assert.That(MslCrc32.Compute(data), Is.EqualTo(expected));

		uint crc = MslCrc32.Initial;
		int position = 0;
		while (position < data.Length)
		{
			int piece = Math.Min(data.Length - position, random.Next(1, 5000));
			crc = MslCrc32.Update(crc, data.AsSpan(position, piece));
			position += piece;
		}
		Assert.That(MslCrc32.Finish(crc), Is.EqualTo(expected));
	}

	private static uint BitwiseCrc32(ReadOnlySpan<byte> data)
	{
		uint crc = 0xFFFF_FFFFu;
		foreach (byte b in data)
		{
			crc ^= b;
			for (int bit = 0; bit < 8; bit++)
				crc = (crc & 1u) != 0 ? (crc >> 1) ^ 0xEDB8_8320u : crc >> 1;
		}
		return ~crc;
	}

	// ── Review fixes ──────────────────────────────────────────────────────────

	private const double CustomLoss = -203.0794;  // HexNAc: not a named loss, so stored in the ext table

	private static List<MslLibraryEntry> MakeEntriesWithCustomLoss(int count)
	{
		List<MslLibraryEntry> entries = MakeEntries(count);
		foreach (MslLibraryEntry e in entries)
			e.MatchedFragmentIons[0].NeutralLoss = CustomLoss;
		return entries;
	}

	private static void AssertCustomLossesRead(MslLibraryData data)
	{
		for (int i = 0; i < data.Count; i++)
		{
			List<MslFragmentIon> ions = data.IsIndexOnly ? data.LoadFragmentsOnDemand(i) : data.Entries[i].MatchedFragmentIons;
			Assert.That(ions[0].NeutralLoss, Is.EqualTo(CustomLoss).Within(1e-9), $"entry {i}");
		}
	}

	/// <summary>
	/// The int32 header field cannot hold the custom-loss table's offset in a file over 2 GB, so
	/// the writer stores 0 there and readers locate the table from the layout. Zeroing the field
	/// (what a &gt;2 GB file carries) or filling it with garbage must not change what is read.
	/// </summary>
	[TestCase(0, 0)]
	[TestCase(0, 12345)]
	[TestCase(3, 0)]
	public void CustomLossTable_IsFoundFromLayout_NotFromHeaderField(int compressionLevel, int headerFieldValue)
	{
		string path = TempPath($"{nameof(CustomLossTable_IsFoundFromLayout_NotFromHeaderField)}_{compressionLevel}_{headerFieldValue}");
		MslWriter.Write(path, MakeEntriesWithCustomLoss(40), compressionLevel);

		byte[] bytes = File.ReadAllBytes(path);
		BitConverter.GetBytes(headerFieldValue).CopyTo(bytes, 28);  // ExtAnnotationTableOffset
		RecomputeCrc(bytes);
		File.WriteAllBytes(path, bytes);

		AssertCustomLossesRead(MslReader.Load(path));
		using MslLibraryData indexOnly = MslReader.LoadIndexOnly(path);
		AssertCustomLossesRead(indexOnly);
	}

	/// <summary>
	/// The located table must end exactly where the offset table begins; a corrupt count is
	/// rejected rather than read as masses.
	/// </summary>
	[Test]
	public void CustomLossTable_CountDisagreesWithLayout_Throws()
	{
		string path = TempPath(nameof(CustomLossTable_CountDisagreesWithLayout_Throws));
		MslWriter.Write(path, MakeEntriesWithCustomLoss(10));

		byte[] bytes = File.ReadAllBytes(path);
		int tableStart = BitConverter.ToInt32(bytes, 28);
		int count = BitConverter.ToInt32(bytes, tableStart);
		BitConverter.GetBytes(count + 1).CopyTo(bytes, tableStart);
		RecomputeCrc(bytes);
		File.WriteAllBytes(path, bytes);

		Assert.That(() => MslReader.Load(path), Throws.TypeOf<FormatException>());
		Assert.That(() => MslReader.LoadIndexOnly(path), Throws.TypeOf<FormatException>());
	}

	/// <summary>
	/// When writing compressed fragments fails part-way, the caller must see the original
	/// exception, not the zstd "pledged size" error from closing the incomplete frame, and no
	/// temp files may be left behind.
	/// </summary>
	[Test]
	public void CompressedWrite_FailurePartWay_SurfacesOriginalException()
	{
		string path = TempPath(nameof(CompressedWrite_FailurePartWay_SurfacesOriginalException));
		List<MslLibraryEntry> entries = MakeEntries(20);
		entries[10].MatchedFragmentIons = null!;  // the layout tolerates null; writing records does not

		Assert.That(() => MslWriter.Write(path, entries, compressionLevel: 3), Throws.TypeOf<NullReferenceException>());
		Assert.That(Directory.GetFiles(OutputDirectory, Path.GetFileName(path) + "*"), Is.Empty);
	}

	/// <summary>
	/// TotalBodyBytes in the string-table header is informational: strings are parsed from their
	/// length prefixes, so a writer that leaves it 0 still produces a readable file in both modes.
	/// </summary>
	[Test]
	public void StringTable_TotalBodyBytesZero_StillReads()
	{
		string path = TempPath(nameof(StringTable_TotalBodyBytesZero_StillReads));
		List<MslLibraryEntry> written = MakeEntries(30);
		MslWriter.Write(path, written);

		long stringTable = MslReader.ReadHeaderOnly(path).StringTableOffset;
		byte[] bytes = File.ReadAllBytes(path);
		BitConverter.GetBytes(0).CopyTo(bytes, (int)stringTable + 4);
		RecomputeCrc(bytes);
		File.WriteAllBytes(path, bytes);

		AssertSameFragments(written, MslReader.Load(path).Entries);
		using MslLibraryData indexOnly = MslReader.LoadIndexOnly(path);
		Assert.That(indexOnly.Entries.Select(e => e.FullSequence), Is.EquivalentTo(written.Select(e => e.FullSequence)));
	}

	/// <summary>
	/// The compressed section must decode to exactly the uncompressed size in the descriptor.
	/// Declaring less or more than the frame holds is rejected, as the old whole-buffer decode did.
	/// </summary>
	[TestCase(-20)]
	[TestCase(20)]
	public void CompressedLoad_DeclaredUncompressedSizeWrong_Throws(int delta)
	{
		string path = TempPath($"{nameof(CompressedLoad_DeclaredUncompressedSizeWrong_Throws)}_{delta}");
		MslWriter.Write(path, MakeEntries(50), compressionLevel: 3);

		byte[] bytes = File.ReadAllBytes(path);
		long declared = BitConverter.ToInt64(bytes, MslFormat.HeaderSize + 8);
		BitConverter.GetBytes(declared + delta).CopyTo(bytes, MslFormat.HeaderSize + 8);
		RecomputeCrc(bytes);
		File.WriteAllBytes(path, bytes);

		Assert.That(() => MslReader.Load(path), Throws.Exception);
	}

	/// <summary>
	/// Files that are not memory-mapped (network or removable drives) are read with positional
	/// reads; that path must match the full load under concurrency, and fail cleanly after Dispose.
	/// </summary>
	[Test]
	public void LoadIndexOnly_PositionalReads_ParallelMatchFullLoad()
	{
		string path = TempPath(nameof(LoadIndexOnly_PositionalReads_ParallelMatchFullLoad));
		List<MslLibraryEntry> written = MakeEntriesWithCustomLoss(300);
		MslWriter.Write(path, written);

		MslLibraryData full = MslReader.Load(path);
		MslLibraryData positional = MslReader.LoadIndexOnly(path, memoryMapFragments: false);

		int mismatches = 0;
		Parallel.For(0, 4 * written.Count, new ParallelOptions { MaxDegreeOfParallelism = 8 }, k =>
		{
			int i = k % written.Count;
			if (!positional.LoadFragmentsOnDemand(i).Select(f => (float)f.Mz).SequenceEqual(FragmentMzs(full.Entries[i])))
				Interlocked.Increment(ref mismatches);
		});
		Assert.That(mismatches, Is.EqualTo(0));
		AssertCustomLossesRead(positional);

		positional.Dispose();
		Assert.That(() => positional.LoadFragmentsOnDemand(0), Throws.TypeOf<ObjectDisposedException>());
	}

	/// <summary>
	/// Mapped reads hold no lock, so Dispose can run while other threads read. Every read must
	/// then either return the right fragments or throw <see cref="ObjectDisposedException"/>;
	/// nothing may touch the unmapped view.
	/// </summary>
	[Test]
	public void LoadIndexOnly_MappedReads_DisposeDuringParallelReads_ReadOrThrowDisposed()
	{
		string path = TempPath(nameof(LoadIndexOnly_MappedReads_DisposeDuringParallelReads_ReadOrThrowDisposed));
		List<MslLibraryEntry> written = MakeEntries(300);
		MslWriter.Write(path, written);
		float[][] expected = MslReader.Load(path).Entries.Select(FragmentMzs).ToArray();

		MslLibraryData mapped = MslReader.LoadIndexOnly(path, memoryMapFragments: true);
		Assert.That(mapped.IsMemoryMapped, Is.True);

		const int workers = 8;
		const int maxReadsPerWorker = 5_000_000;
		long reads = 0;
		int mismatches = 0, otherErrors = 0, sawDisposed = 0;
		Task[] tasks = Enumerable.Range(0, workers).Select(w => Task.Run(() =>
		{
			for (int k = 0; k < maxReadsPerWorker; k++)
			{
				int i = (w * 37 + k) % expected.Length;
				try
				{
					if (!mapped.LoadFragmentsOnDemand(i).Select(f => (float)f.Mz).SequenceEqual(expected[i]))
						Interlocked.Increment(ref mismatches);
					Interlocked.Increment(ref reads);
				}
				catch (ObjectDisposedException)
				{
					Interlocked.Increment(ref sawDisposed);
					return;
				}
				catch
				{
					Interlocked.Increment(ref otherErrors);
					return;
				}
			}
		})).ToArray();

		SpinWait.SpinUntil(() => Interlocked.Read(ref reads) >= 10_000, TimeSpan.FromSeconds(30));
		mapped.Dispose();
		Assert.That(Task.WaitAll(tasks, TimeSpan.FromSeconds(60)), Is.True, "readers did not stop");

		Assert.That(Interlocked.Read(ref reads), Is.GreaterThanOrEqualTo(10_000), "reads had not started");
		Assert.That(mismatches, Is.EqualTo(0));
		Assert.That(otherErrors, Is.EqualTo(0));
		Assert.That(sawDisposed, Is.EqualTo(workers), "every reader ends on ObjectDisposedException");
	}

	/// <summary>
	/// When the OS refuses the memory map, the index-only load falls back to positional reads
	/// instead of failing.
	/// </summary>
	[TestCase(typeof(IOException))]
	[TestCase(typeof(UnauthorizedAccessException))]
	public void LoadIndexOnly_MapCannotBeCreated_FallsBackToPositionalReads(Type exceptionType)
	{
		string path = TempPath(nameof(LoadIndexOnly_MapCannotBeCreated_FallsBackToPositionalReads) + exceptionType.Name);
		List<MslLibraryEntry> written = MakeEntriesWithCustomLoss(50);
		MslWriter.Write(path, written);
		MslLibraryData full = MslReader.Load(path);

		MslLibraryData.CreateMapForTesting = _ => throw (Exception)Activator.CreateInstance(exceptionType, "no map")!;
		MslLibraryData fallback;
		try
		{
			fallback = MslReader.LoadIndexOnly(path, memoryMapFragments: true);
		}
		finally
		{
			MslLibraryData.CreateMapForTesting = null;
		}

		using (fallback)
		{
			Assert.That(fallback.IsMemoryMapped, Is.False);
			for (int i = 0; i < written.Count; i++)
				Assert.That(fallback.LoadFragmentsOnDemand(i).Select(f => (float)f.Mz), Is.EqualTo(FragmentMzs(full.Entries[i])));
			AssertCustomLossesRead(fallback);
		}
	}

	/// <summary>
	/// Any other failure while mapping still fails the load.
	/// </summary>
	[Test]
	public void LoadIndexOnly_MapFailsWithUnexpectedException_Throws()
	{
		string path = TempPath(nameof(LoadIndexOnly_MapFailsWithUnexpectedException_Throws));
		MslWriter.Write(path, MakeEntries(5));

		MslLibraryData.CreateMapForTesting = _ => throw new InvalidOperationException("unexpected");
		try
		{
			Assert.That(() => MslReader.LoadIndexOnly(path, memoryMapFragments: true),
				Throws.TypeOf<InvalidOperationException>());
		}
		finally
		{
			MslLibraryData.CreateMapForTesting = null;
		}
		Assert.That(() => File.Delete(path), Throws.Nothing, "the failed load released the file");
	}

	/// <summary>
	/// The precursor section is read in chunks so the byte span never overflows; a chunk size
	/// that does not divide the count must give the same records as one read.
	/// </summary>
	[Test]
	public void ReadPrecursorArray_InChunks_MatchesSingleRead()
	{
		string path = TempPath(nameof(ReadPrecursorArray_InChunks_MatchesSingleRead));
		MslWriter.Write(path, MakeEntries(300));
		MslFileHeader header;
		using (MslLibraryData data = MslReader.LoadIndexOnly(path))
			header = data.Header;

		using FileStream fs = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.Read);
		MslPrecursorRecord[] whole = MslReader.ReadPrecursorArrayFromStream(fs, header, chunkRecords: int.MaxValue);
		MslPrecursorRecord[] chunked = MslReader.ReadPrecursorArrayFromStream(fs, header, chunkRecords: 7);

		Assert.That(chunked.Length, Is.EqualTo(300));
		Assert.That(System.Runtime.InteropServices.MemoryMarshal.AsBytes(chunked.AsSpan()).SequenceEqual(
			System.Runtime.InteropServices.MemoryMarshal.AsBytes(whole.AsSpan())), Is.True);
	}

	/// <summary>
	/// A UNC path is a network share, so it is never memory-mapped; an ordinary path is decided
	/// without throwing.
	/// </summary>
	[Test]
	public void IsOnLocalFixedDrive_UncPathIsNotLocal()
	{
		Assert.That(MslReader.IsOnLocalFixedDrive(@"\\server\share\library.msl"), Is.False);
		Assert.That(() => MslReader.IsOnLocalFixedDrive(TempPath("any")), Throws.Nothing);
	}

	private static void RecomputeCrc(byte[] bytes)
	{
		int footerStart = bytes.Length - MslFormat.FooterSize;
		long offsetTableOffset = BitConverter.ToInt64(bytes, footerStart);
		uint newCrc = MslWriter.ComputeCrc32OfArray(bytes, (int)offsetTableOffset);
		BitConverter.GetBytes(newCrc).CopyTo(bytes, footerStart + 12);
	}
}
