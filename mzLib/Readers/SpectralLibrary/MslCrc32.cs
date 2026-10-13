using System.Buffers.Binary;

namespace Readers.SpectralLibrary;

/// <summary>
/// CRC-32/ISO-HDLC (reflected polynomial 0xEDB88320; the same result as zlib <c>crc32()</c>
/// and PKZIP), shared by <see cref="MslWriter"/> and <see cref="MslReader"/>.
/// <para>
/// Computed eight bytes per step with slicing-by-8 tables. The result is identical to the
/// one-byte-at-a-time table loop, several times faster; the CRC pass over the data section
/// is the largest single cost of opening a library.
/// </para>
/// </summary>
internal static class MslCrc32
{
	/// <summary>Initial value of the CRC register.</summary>
	internal const uint Initial = 0xFFFF_FFFFu;

	/// <summary>
	/// Eight 256-entry tables, back to back. Table 0 is the classic byte table; table k
	/// advances a byte by k further zero bytes.
	/// </summary>
	private static readonly uint[] Tables = BuildTables();

	private static uint[] BuildTables()
	{
		const uint Polynomial = 0xEDB8_8320u;
		var t = new uint[8 * 256];

		for (uint i = 0; i < 256; i++)
		{
			uint entry = i;
			for (int bit = 0; bit < 8; bit++)
				entry = (entry & 1u) != 0 ? (entry >> 1) ^ Polynomial : entry >> 1;
			t[i] = entry;
		}

		for (int k = 1; k < 8; k++)
			for (int i = 0; i < 256; i++)
			{
				uint previous = t[(k - 1) * 256 + i];
				t[k * 256 + i] = (previous >> 8) ^ t[previous & 0xFF];
			}

		return t;
	}

	/// <summary>
	/// Feeds <paramref name="data"/> into a running CRC register. Start from
	/// <see cref="Initial"/> and finish with <see cref="Finish"/>.
	/// </summary>
	internal static uint Update(uint crc, ReadOnlySpan<byte> data)
	{
		ReadOnlySpan<uint> t = Tables;
		int i = 0;

		for (; i + 8 <= data.Length; i += 8)
		{
			uint one = BinaryPrimitives.ReadUInt32LittleEndian(data.Slice(i)) ^ crc;
			uint two = BinaryPrimitives.ReadUInt32LittleEndian(data.Slice(i + 4));
			crc = t[7 * 256 + (int)(one & 0xFF)]
				^ t[6 * 256 + (int)((one >> 8) & 0xFF)]
				^ t[5 * 256 + (int)((one >> 16) & 0xFF)]
				^ t[4 * 256 + (int)(one >> 24)]
				^ t[3 * 256 + (int)(two & 0xFF)]
				^ t[2 * 256 + (int)((two >> 8) & 0xFF)]
				^ t[1 * 256 + (int)((two >> 16) & 0xFF)]
				^ t[(int)(two >> 24)];
		}

		for (; i < data.Length; i++)
			crc = (crc >> 8) ^ t[(int)((crc ^ data[i]) & 0xFF)];

		return crc;
	}

	/// <summary>Final one's complement that turns the register into the checksum.</summary>
	internal static uint Finish(uint crc) => crc ^ 0xFFFF_FFFFu;

	/// <summary>Checksum of a whole span.</summary>
	internal static uint Compute(ReadOnlySpan<byte> data) => Finish(Update(Initial, data));
}
