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

	private static void RecomputeCrc(byte[] bytes)
	{
		int footerStart = bytes.Length - MslFormat.FooterSize;
		long offsetTableOffset = BitConverter.ToInt64(bytes, footerStart);
		uint newCrc = MslWriter.ComputeCrc32OfArray(bytes, (int)offsetTableOffset);
		BitConverter.GetBytes(newCrc).CopyTo(bytes, footerStart + 12);
	}
}
