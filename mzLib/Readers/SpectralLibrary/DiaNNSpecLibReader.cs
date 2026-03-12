using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using Omics.SpectrumMatch;

namespace Readers.SpectralLibrary.DiaNNSpectralLibrary
{
	/// <summary>
	/// Reads DIA-NN 2.3.2 .speclib binary spectral library files (format version −10).
	///
	/// The format is anchor-based: there is no fixed-size precursor record array.
	/// Every precursor is located by scanning the raw file bytes for one of six known
	/// 4-byte marker patterns. The binary specification is fully documented in
	/// Session7_SpecLibWrapup.md.
	///
	/// Reading strategy:
	///   1. Load the entire file into a byte array (avoids repeated seeking).
	///   2. Parse the header and protein table sequentially.
	///   3. Scan all bytes from end-of-protein-table to EOF for all six marker patterns.
	///   4. For each marker hit, validate the surrounding context and extract the record.
	///   5. Sort and return results.
	///
	/// Fragment annotation (ion series, ion number) is intentionally NOT decoded.
	/// Bytes 0–2 of each fragment record are opaque DIA-NN internal indices that are
	/// not reliable indicators of ion type. Only fragment_mz and fragment_intensity
	/// are extracted.
	/// </summary>
	public static class DiaNNSpecLibReader
	{
		// ─── UniMod tag pattern ──────────────────────────────────────────────────────
		private static readonly Regex UniModTagRegex =
			new Regex(@"\(UniMod:\d+\)", RegexOptions.Compiled);

		// ─── Public API ──────────────────────────────────────────────────────────────

		/// <summary>
		/// Reads all precursor entries from a DIA-NN .speclib binary file.
		/// </summary>
		/// <param name="filePath">Path to the .speclib file.</param>
		/// <returns>
		/// List of <see cref="DiaNNLibraryEntry"/> objects, sorted by global_index,
		/// each populated with protein metadata, precursor fields, and fragment ions.
		/// Fragment null-placeholder slots (mz = 0.0) are discarded.
		/// </returns>
		/// <exception cref="FileNotFoundException">File does not exist.</exception>
		/// <exception cref="FormatException">File is not a valid DIA-NN 2.3.2 .speclib.</exception>
		public static List<DiaNNLibraryEntry> Read(string filePath)
		{
			if (!File.Exists(filePath))
				throw new FileNotFoundException($"DIA-NN .speclib file not found: {filePath}", filePath);

			byte[] data = File.ReadAllBytes(filePath);

			int pos = 0;

			// ── 1. Header ─────────────────────────────────────────────────────────────
			int version = ReadInt32(data, pos); pos += 4;
			if (version != DiaNNBinaryStructs.ExpectedVersion)
				throw new FormatException(
					$"Unexpected .speclib version {version}. Expected {DiaNNBinaryStructs.ExpectedVersion}. " +
					$"File: {filePath}");

			// Four opaque flag int32 fields
			pos += 4 * 4;

			// FASTA path (length-prefixed, UTF-8, NOT null-terminated)
			int fastaPathLen = ReadInt32(data, pos); pos += 4;
			string fastaPath = Encoding.UTF8.GetString(data, pos, fastaPathLen);
			pos += fastaPathLen;

			// ── 2. Protein table ──────────────────────────────────────────────────────
			int nProteins = ReadInt32(data, pos); pos += 4;

			var proteins = new List<ProteinInfo>(nProteins);
			for (int i = 0; i < nProteins; i++)
			{
				// Protein 0: 2 opaque int32 prefix fields
				// Proteins 1..N-1: 4 opaque int32 prefix fields
				int prefixCount = (i == 0) ? 2 : 4;
				pos += prefixCount * 4;

				string accession = ReadLengthPrefixedString(data, ref pos);
				string protName = ReadLengthPrefixedString(data, ref pos);
				string geneName = ReadLengthPrefixedString(data, ref pos);

				proteins.Add(new ProteinInfo(i, accession, protName, geneName));
			}

			// ── 3. Scan for precursor markers ─────────────────────────────────────────
			// Skip the secondary table — scan from the current position forward.
			// The scanner validates each hit before accepting it, so false positives
			// in the secondary table are harmless.
			var entries = ScanForPrecursors(data, pos, proteins);

			entries.Sort((a, b) => a.GlobalIndex.CompareTo(b.GlobalIndex));

			return entries.Select(e => e.ToLibraryEntry()).ToList();
		}

		// ─── Marker scanner ───────────────────────────────────────────────────────────

		private static List<PrecursorRecord> ScanForPrecursors(
			byte[] data, int searchStart, List<ProteinInfo> proteins)
		{
			var results = new List<PrecursorRecord>(512);
			int dataLen = data.Length;

			// We need at least 48 bytes before the marker (pre-name block)
			// and enough bytes after for a minimal name + post-name block.
			// marker_pos must satisfy: marker_pos >= 48, marker_pos + 8 + 4 + 5 < dataLen
			int minPos = Math.Max(searchStart, 48);
			int maxPos = dataLen - 4;

			for (int i = minPos; i < maxPos - 2; i++)
			{
				// Fast pre-check: all markers start with 0x07, 0x00
				if (data[i] != 0x07 || data[i + 1] != 0x00)
					continue;

				var terminusType = DiaNNBinaryStructs.ClassifyMarker(
					data[i], data[i + 1], data[i + 2], data[i + 3]);
				if (terminusType == null)
					continue;

				int markerPos = i;
				int nlp = markerPos + 8;  // position of name_length int32

				if (nlp + 4 >= dataLen)
					continue;

				// Read pre-name fields
				int globalIndex = ReadInt32(data, nlp + DiaNNBinaryStructs.OffsetGlobalIndex);
				int charge = ReadInt32(data, nlp + DiaNNBinaryStructs.OffsetCharge);
				int strippedSeqLen = ReadInt32(data, nlp + DiaNNBinaryStructs.OffsetStrippedSeqLen);
				float precursorMz = ReadFloat(data, nlp + DiaNNBinaryStructs.OffsetPrecursorMz);
				float iRT = ReadFloat(data, nlp + DiaNNBinaryStructs.OffsetIRT);
				float ionMobility = ReadFloat(data, nlp + DiaNNBinaryStructs.OffsetIonMobility);
				int protGroupIndex = ReadInt32(data, nlp + DiaNNBinaryStructs.OffsetProteinGroupIndex);

				// Validation gate
				if (globalIndex < 0 || globalIndex >= DiaNNBinaryStructs.GlobalIndexMax)
					continue;
				if (precursorMz < DiaNNBinaryStructs.PrecursorMzMin ||
					precursorMz > DiaNNBinaryStructs.PrecursorMzMax)
					continue;
				if (charge < 1 || charge > 4)
					continue;

				// Name field
				int nameLength = ReadInt32(data, nlp);
				if (nameLength < DiaNNBinaryStructs.NameLengthMin ||
					nameLength > DiaNNBinaryStructs.NameLengthMax)
					continue;

				int nameStart = nlp + 4;
				if (nameStart + nameLength > dataLen)
					continue;

				// Validate name bytes: all printable ASCII, last byte is digit '1'–'4'
				if (!IsValidNameBytes(data, nameStart, nameLength))
					continue;

				string name = Encoding.ASCII.GetString(data, nameStart, nameLength);

				// Post-name block
				int nameEnd = nameStart + nameLength;
				if (nameEnd + DiaNNBinaryStructs.PostNameBlockSize > dataLen)
					continue;

				int nFragments = ReadInt32(data, nameEnd + DiaNNBinaryStructs.OffsetNFragments);
				if (nFragments < DiaNNBinaryStructs.NFragmentsMin ||
					nFragments > DiaNNBinaryStructs.NFragmentsMax)
					continue;

				float topFragMz = ReadFloat(data, nameEnd + DiaNNBinaryStructs.OffsetTopFragmentMz);
				float topFragInty = ReadFloat(data, nameEnd + DiaNNBinaryStructs.OffsetTopFragmentIntensity);

				// Fragment records
				int fragStart = nameEnd + DiaNNBinaryStructs.PostNameBlockSize;
				int fragBlockSize = nFragments * DiaNNBinaryStructs.FragmentRecordSize;
				if (fragStart + fragBlockSize > dataLen)
					continue;

				var fragments = ReadFragments(data, fragStart, nFragments);

				// Resolve protein from group index (best effort — fall back to first protein)
				ProteinInfo protein = ResolveProtein(proteins, protGroupIndex);

				results.Add(new PrecursorRecord(
					globalIndex, charge, strippedSeqLen,
					precursorMz, iRT, ionMobility,
					terminusType.Value, protGroupIndex,
					name, fragments,
					protein));
			}

			return results;
		}

		// ─── Fragment reading ─────────────────────────────────────────────────────────

		private static List<FragmentData> ReadFragments(byte[] data, int offset, int count)
		{
			var list = new List<FragmentData>(count);
			for (int i = 0; i < count; i++)
			{
				int recOffset = offset + i * DiaNNBinaryStructs.FragmentRecordSize;
				// Bytes 0–3 are opaque — skip.
				float mz = ReadFloat(data, recOffset + DiaNNBinaryStructs.FragmentOffsetMz);
				float intensity = ReadFloat(data, recOffset + DiaNNBinaryStructs.FragmentOffsetIntensity);

				// Discard null placeholder slots
				if (mz == 0.0f)
					continue;

				list.Add(new FragmentData(mz, intensity));
			}
			return list;
		}

		// ─── Name parsing ─────────────────────────────────────────────────────────────

		/// <summary>
		/// Parses a DIA-NN .speclib name string of the form ModifiedSequence + ChargeDigit.
		///
		/// Examples:
		///   "ALGVGLATR2"                  → seq=ALGVGLATR, charge=2, no mods
		///   "AVM(UniMod:35)K1"            → seq=AVMK, charge=1, mods=[(2, 35)]
		///   "(UniMod:1)ASNQTYK2"          → seq=ASNQTYK, charge=2, mods=[(-1, 1)]
		///   "AS(UniMod:21)NQTYK2"         → seq=ASNQTYK, charge=2, mods=[(1, 21)]
		///   "ADVTLAK(UniMod:121)2"        → seq=ADVTLAK, charge=2, mods=[(6, 121)]
		///
		/// Returns the charge digit as int and a list of (position, uniModId) tuples
		/// where position −1 means N-terminal modification.
		/// </summary>
		internal static (string strippedSeq, int charge, List<(int pos, int uniModId)> mods)
			ParseName(string name)
		{
			if (string.IsNullOrEmpty(name) || name.Length < 2)
				throw new FormatException($"Name too short to parse: '{name}'");

			// Last character is the charge digit
			int charge = name[name.Length - 1] - '0';
			string modSeq = name.Substring(0, name.Length - 1);

			var mods = new List<(int pos, int uniModId)>();
			var strippedSb = new StringBuilder(modSeq.Length);

			int idx = 0;
			int aaIndex = 0;  // 0-based index in the stripped sequence

			while (idx < modSeq.Length)
			{
				if (modSeq[idx] == '(')
				{
					// Find closing paren
					int close = modSeq.IndexOf(')', idx);
					if (close < 0)
						throw new FormatException($"Unclosed '(' in name: '{name}'");

					string tag = modSeq.Substring(idx, close - idx + 1);
					int uniModId = ParseUniModId(tag);

					// N-terminal mod: tag appears before any amino acid (aaIndex still 0).
					// Mid/C-terminal mod: tag appears immediately after the residue it modifies,
					// so the position is aaIndex - 1 (the residue just appended).
					int position = (aaIndex == 0) ? -1 : aaIndex - 1;
					mods.Add((position, uniModId));

					idx = close + 1;
				}
				else
				{
					strippedSb.Append(modSeq[idx]);
					aaIndex++;
					idx++;
				}
			}

			return (strippedSb.ToString(), charge, mods);
		}

		private static int ParseUniModId(string tag)
		{
			// tag is "(UniMod:N)"
			const string prefix = "(UniMod:";
			if (!tag.StartsWith(prefix))
				throw new FormatException($"Unexpected modification tag format: '{tag}'");
			string numStr = tag.Substring(prefix.Length, tag.Length - prefix.Length - 1);
			return int.Parse(numStr);
		}

		// ─── Helpers ──────────────────────────────────────────────────────────────────

		private static bool IsValidNameBytes(byte[] data, int start, int length)
		{
			for (int i = 0; i < length; i++)
			{
				byte b = data[start + i];
				if (b < 32 || b > 126) return false;
			}
			// Last byte must be '1'–'4'
			byte last = data[start + length - 1];
			return last >= (byte)'1' && last <= (byte)'4';
		}

		private static int ReadInt32(byte[] data, int offset)
			=> BitConverter.ToInt32(data, offset);

		private static float ReadFloat(byte[] data, int offset)
			=> BitConverter.ToSingle(data, offset);

		private static string ReadLengthPrefixedString(byte[] data, ref int pos)
		{
			int len = ReadInt32(data, pos); pos += 4;
			string s = len > 0 ? Encoding.UTF8.GetString(data, pos, len) : string.Empty;
			pos += len;
			return s;
		}

		private static ProteinInfo ResolveProtein(List<ProteinInfo> proteins, int protGroupIndex)
		{
			// The group index points into the secondary table group list, which we don't parse.
			// For individual groups: DIA-NN places joint groups first, then one group per protein.
			// A safe heuristic: clamp to valid protein index range and use the protein at that index.
			if (proteins.Count == 0) return ProteinInfo.Empty;
			int idx = protGroupIndex < proteins.Count ? protGroupIndex : proteins.Count - 1;
			if (idx < 0) idx = 0;
			return proteins[idx];
		}

		// ─── Internal data models ─────────────────────────────────────────────────────

		private readonly struct FragmentData
		{
			public readonly float Mz;
			public readonly float Intensity;
			public FragmentData(float mz, float intensity) { Mz = mz; Intensity = intensity; }
		}

		private class ProteinInfo
		{
			public readonly int Index;
			public readonly string Accession;
			public readonly string Name;
			public readonly string Gene;

			public ProteinInfo(int index, string accession, string name, string gene)
			{
				Index = index; Accession = accession; Name = name; Gene = gene;
			}

			public static readonly ProteinInfo Empty = new ProteinInfo(0, "", "", "");
		}

		private class PrecursorRecord
		{
			public readonly int GlobalIndex;
			public readonly int Charge;
			public readonly int StrippedSeqLen;
			public readonly float PrecursorMz;
			public readonly float IRT;
			public readonly float IonMobility;
			public readonly DiaNNBinaryStructs.TerminusType TerminusType;
			public readonly int ProteinGroupIndex;
			public readonly string Name;
			public readonly List<FragmentData> Fragments;
			public readonly ProteinInfo Protein;

			public PrecursorRecord(
				int globalIndex, int charge, int strippedSeqLen,
				float precursorMz, float irt, float ionMobility,
				DiaNNBinaryStructs.TerminusType terminusType, int protGroupIndex,
				string name, List<FragmentData> fragments, ProteinInfo protein)
			{
				GlobalIndex = globalIndex;
				Charge = charge;
				StrippedSeqLen = strippedSeqLen;
				PrecursorMz = precursorMz;
				IRT = irt;
				IonMobility = ionMobility;
				TerminusType = terminusType;
				ProteinGroupIndex = protGroupIndex;
				Name = name;
				Fragments = fragments;
				Protein = protein;
			}

			public DiaNNLibraryEntry ToLibraryEntry()
			{
				var (strippedSeq, charge, mods) = ParseName(Name);

				// Build DIA-NN format modified sequence (parenthesis notation, no underscores)
				string modifiedSequence = BuildModifiedSequence(strippedSeq, mods);

				var fragmentIons = Fragments.Select(f => new DiaNNFragmentIon
				{
					Mz = f.Mz,
					Intensity = f.Intensity,
					// Ion type annotation is not decoded (bytes 0-2 are opaque)
					IonType = 'y',           // placeholder
					SeriesNumber = 0,        // placeholder
					Charge = 1,              // placeholder
					LossType = "noloss",
				}).ToList();

				return new DiaNNLibraryEntry
				{
					ModifiedSequence = modifiedSequence,
					StrippedSequence = strippedSeq,
					PrecursorMz = PrecursorMz,
					PrecursorCharge = Charge,
					RetentionTime = IRT,
					IonMobility = IonMobility,
					IsDecoy = false,
					Fragments = fragmentIons,
					ProteinId = Protein.Accession,
					ProteinName = Protein.Name,
					GeneName = Protein.Gene,
					IsProteotypic = TerminusType != DiaNNBinaryStructs.TerminusType.SharedConflict
									 && TerminusType != DiaNNBinaryStructs.TerminusType.SharedSameType,
					TerminusType = TerminusType,
					ProteinGroupIndex = ProteinGroupIndex,
					GlobalIndex = GlobalIndex,
				};
			}

			/// <summary>
			/// Reconstructs the DIA-NN parenthesis-notation modified sequence from
			/// a stripped sequence and a list of (position, uniModId) tuples.
			/// N-terminal modifications have position = −1.
			/// </summary>
			private static string BuildModifiedSequence(
				string stripped, List<(int pos, int uniModId)> mods)
			{
				var sb = new StringBuilder(stripped.Length + mods.Count * 12);

				// N-terminal mod (position = -1)
				foreach (var (pos, uid) in mods.Where(m => m.pos == -1))
					sb.Append($"(UniMod:{uid})");

				for (int i = 0; i < stripped.Length; i++)
				{
					foreach (var (pos, uid) in mods.Where(m => m.pos == i))
						sb.Append($"(UniMod:{uid})");
					sb.Append(stripped[i]);
				}

				return sb.ToString();
			}
		}
	}
}