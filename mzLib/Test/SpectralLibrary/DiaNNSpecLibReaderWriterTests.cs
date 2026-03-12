using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using NUnit.Framework;

namespace Readers.SpectralLibrary.DiaNNSpectralLibrary.Tests
{
	/// <summary>
	/// Critical unit tests for the rewritten DIA-NN .speclib binary reader and writer.
	///
	/// Each test documents:
	///   WHAT: the specific behaviour under test.
	///   WHY:  the consequence if this behaviour is wrong.
	///
	/// Tests are deliberately minimal. We cover only the behaviours where a bug would
	/// silently produce wrong results (silent data loss, wrong field values, format
	/// corruption). Edge-case and boundary tests are omitted unless they correspond to
	/// a known pitfall documented in the session findings.
	/// </summary>
	[TestFixture]
	public class DiaNNSpecLibReaderWriterTests
	{
		// ═══════════════════════════════════════════════════════════════════════════════
		// PART 1 — DiaNNBinaryStructs: marker classification
		// ═══════════════════════════════════════════════════════════════════════════════

		/// <summary>
		/// WHAT: ClassifyMarker correctly identifies all six marker byte patterns.
		/// WHY:  A scanner that misses even one pattern silently drops entire categories
		///       of precursors (e.g. all shared-peptide entries if the subnormal patterns
		///       are not recognised). Verified against Session4/Session7 findings.
		/// </summary>
		[Test]
		public void ClassifyMarker_AllSixPatterns_ReturnCorrectTerminusType()
		{
			Assert.That(DiaNNBinaryStructs.ClassifyMarker(0x07, 0x00, 0x3C, 0x40),
				Is.EqualTo(DiaNNBinaryStructs.TerminusType.Internal));

			Assert.That(DiaNNBinaryStructs.ClassifyMarker(0x07, 0x00, 0x39, 0x40),
				Is.EqualTo(DiaNNBinaryStructs.TerminusType.NTerminal));

			Assert.That(DiaNNBinaryStructs.ClassifyMarker(0x07, 0x00, 0x36, 0x40),
				Is.EqualTo(DiaNNBinaryStructs.TerminusType.CTerminal));

			Assert.That(DiaNNBinaryStructs.ClassifyMarker(0x07, 0x00, 0x33, 0x40),
				Is.EqualTo(DiaNNBinaryStructs.TerminusType.NAndCTerminal));

			// Subnormal floats — the critical shared-peptide markers that old scanners miss.
			Assert.That(DiaNNBinaryStructs.ClassifyMarker(0x07, 0x00, 0x3D, 0x00),
				Is.EqualTo(DiaNNBinaryStructs.TerminusType.SharedConflict));

			Assert.That(DiaNNBinaryStructs.ClassifyMarker(0x07, 0x00, 0x39, 0x00),
				Is.EqualTo(DiaNNBinaryStructs.TerminusType.SharedSameType));
		}

		/// <summary>
		/// WHAT: ClassifyMarker returns null for byte patterns that are not markers.
		/// WHY:  Random data bytes will inevitably start with 0x07, 0x00. The scanner
		///       relies on a null return to skip false positives, so incorrect acceptance
		///       would fabricate phantom precursors.
		/// </summary>
		[Test]
		public void ClassifyMarker_UnknownPattern_ReturnsNull()
		{
			Assert.That(DiaNNBinaryStructs.ClassifyMarker(0x07, 0x00, 0x3C, 0x00), Is.Null,
				"Pattern similar to SharedConflict but wrong third byte should not match.");

			Assert.That(DiaNNBinaryStructs.ClassifyMarker(0x07, 0x01, 0x3C, 0x40), Is.Null,
				"Wrong second byte should not match.");

			Assert.That(DiaNNBinaryStructs.ClassifyMarker(0xFF, 0x00, 0x3C, 0x40), Is.Null,
				"Wrong first byte should not match.");
		}

		// ═══════════════════════════════════════════════════════════════════════════════
		// PART 2 — ParseName: name string parsing
		// ═══════════════════════════════════════════════════════════════════════════════

		/// <summary>
		/// WHAT: Unmodified peptide name parses to correct stripped sequence, charge, and empty mod list.
		/// WHY:  ParseName is the primary way sequences are extracted from the binary file.
		///       A bug here corrupts every sequence in the output library.
		/// </summary>
		[Test]
		public void ParseName_UnmodifiedPeptide_ReturnsCorrectFields()
		{
			var (seq, charge, mods) = DiaNNSpecLibReader.ParseName("ALGVGLATR2");

			Assert.That(seq, Is.EqualTo("ALGVGLATR"));
			Assert.That(charge, Is.EqualTo(2));
			Assert.That(mods, Is.Empty);
		}

		/// <summary>
		/// WHAT: A mid-sequence modification is assigned the correct 0-based position in
		///       the stripped sequence.
		/// WHY:  Positions are used to reconstruct the modified sequence and to map
		///       modifications onto the peptide for MetaMorpheus scoring. Off-by-one errors
		///       here would silently misplace every modification.
		/// </summary>
		[Test]
		public void ParseName_MidSequenceModification_CorrectPosition()
		{
			// "AVM(UniMod:35)K1" — oxidation on M at position 2 (0-based in AVMK)
			var (seq, charge, mods) = DiaNNSpecLibReader.ParseName("AVM(UniMod:35)K1");

			Assert.That(seq, Is.EqualTo("AVMK"));
			Assert.That(charge, Is.EqualTo(1));
			Assert.That(mods, Has.Count.EqualTo(1));
			Assert.That(mods[0].pos, Is.EqualTo(2), "Oxidised M is at 0-based position 2.");
			Assert.That(mods[0].uniModId, Is.EqualTo(35));
		}

		/// <summary>
		/// WHAT: A tag that appears before any amino acid letter is classified as N-terminal
		///       (position = −1), not as a modification on the first residue (position = 0).
		/// WHY:  N-terminal modifications are a distinct concept from residue modifications.
		///       Misclassifying them as position 0 would assign acetylation to the first
		///       amino acid rather than to the peptide N-terminus, breaking downstream
		///       mass calculation.
		/// </summary>
		[Test]
		public void ParseName_NTerminalModification_PositionIsMinusOne()
		{
			// "(UniMod:1)ASNQTYK2" — N-terminal acetylation
			var (seq, charge, mods) = DiaNNSpecLibReader.ParseName("(UniMod:1)ASNQTYK2");

			Assert.That(seq, Is.EqualTo("ASNQTYK"));
			Assert.That(charge, Is.EqualTo(2));
			Assert.That(mods, Has.Count.EqualTo(1));
			Assert.That(mods[0].pos, Is.EqualTo(-1), "N-terminal mod must have position −1.");
			Assert.That(mods[0].uniModId, Is.EqualTo(1));
		}

		/// <summary>
		/// WHAT: A modification on the last residue before the charge digit is parsed correctly.
		/// WHY:  "ADVTLAK(UniMod:121)2" has the tag after the final amino acid K. The parser
		///       must not confuse the tag's closing ')' with the charge digit. GlyGly (121) on
		///       C-terminal K is common in ubiquitin proteomics.
		/// </summary>
		[Test]
		public void ParseName_ModificationOnLastResidue_ParsedCorrectly()
		{
			// GlyGly on K at position 6 (0-based in ADVTLAK, length 7)
			var (seq, charge, mods) = DiaNNSpecLibReader.ParseName("ADVTLAK(UniMod:121)2");

			Assert.That(seq, Is.EqualTo("ADVTLAK"));
			Assert.That(charge, Is.EqualTo(2));
			Assert.That(mods, Has.Count.EqualTo(1));
			Assert.That(mods[0].pos, Is.EqualTo(6));
			Assert.That(mods[0].uniModId, Is.EqualTo(121));
		}

		// ═══════════════════════════════════════════════════════════════════════════════
		// PART 3 — Round-trip: Write then Read produces identical data
		// ═══════════════════════════════════════════════════════════════════════════════

		/// <summary>
		/// Builds a minimal but valid library entry for round-trip tests.
		/// Uses real m/z and iRT values from the synthetic2 MSP file (Session 7).
		/// </summary>
		private static DiaNNLibraryEntry MakeEntry(
			int globalIndex,
			string modSeq,
			string strippedSeq,
			int charge,
			float mz,
			float irt,
			DiaNNBinaryStructs.TerminusType terminus,
			string accession = "TEST_PROT",
			params (float mz, float intensity)[] fragments)
		{
			var frags = fragments.Select(f => new DiaNNFragmentIon
			{
				Mz = f.mz,
				Intensity = f.intensity,
				IonType = 'y',
				SeriesNumber = 1,
				Charge = 1,
				LossType = "noloss",
			}).ToList();

			return new DiaNNLibraryEntry
			{
				GlobalIndex = globalIndex,
				ModifiedSequence = modSeq,
				StrippedSequence = strippedSeq,
				PrecursorCharge = charge,
				PrecursorMz = mz,
				RetentionTime = irt,
				IonMobility = 0.0,
				TerminusType = terminus,
				ProteinGroupIndex = 0,
				ProteinId = accession,
				ProteinName = accession + "_NAME",
				GeneName = "GENE1",
				IsDecoy = false,
				IsProteotypic = true,
				Fragments = frags,
			};
		}

		/// <summary>
		/// WHAT: A library written then read back has the same number of precursors.
		/// WHY:  This is the fundamental correctness invariant. If precursor count changes,
		///       either the writer is not emitting valid marker bytes or the reader's
		///       validation gate is incorrectly rejecting real records.
		/// </summary>
		[Test]
		public void RoundTrip_PrecursorCount_IsPreserved()
		{
			var entries = new List<DiaNNLibraryEntry>
			{
                // Unmodified internal peptide — MSTYNQK/2, real values from Session7 MSP
                MakeEntry(0, "MSTYNQK", "MSTYNQK", 2, 436.2026f, 21.51f,
					DiaNNBinaryStructs.TerminusType.Internal, "PROT_A",
					(740.36f, 1.0f), (653.33f, 0.538f), (552.28f, 0.298f)),

                // N-terminal peptide
                MakeEntry(1, "GILDQSIVR", "GILDQSIVR", 2, 500.793f, 72.48f,
					DiaNNBinaryStructs.TerminusType.NTerminal, "PROT_B",
					(830.47f, 0.98f), (717.39f, 1.0f), (602.36f, 0.557f)),
			};

			string path = Path.Combine(Path.GetTempPath(), $"roundtrip_count_{Guid.NewGuid():N}.speclib");
			try
			{
				DiaNNSpecLibWriter.Write(entries, path);
				var readBack = DiaNNSpecLibReader.Read(path);

				Assert.That(readBack, Has.Count.EqualTo(2),
					"Read-back must return exactly as many precursors as were written.");
			}
			finally
			{
				if (File.Exists(path)) File.Delete(path);
			}
		}

		/// <summary>
		/// WHAT: After a round-trip, precursor m/z and iRT values are preserved within
		///       float precision tolerance.
		/// WHY:  These are the two most critical numerical fields used in DIA window matching
		///       and retention time calibration. Silent corruption here would cause precursors
		///       to be matched against the wrong DIA windows.
		/// </summary>
		[Test]
		public void RoundTrip_PrecursorMzAndIRT_PreservedWithinFloatPrecision()
		{
			float expectedMz = 436.2026f;
			float expectedIRT = 21.51f;

			var entries = new List<DiaNNLibraryEntry>
			{
				MakeEntry(0, "MSTYNQK", "MSTYNQK", 2, expectedMz, expectedIRT,
					DiaNNBinaryStructs.TerminusType.Internal, "PROT_A",
					(740.36f, 1.0f), (653.33f, 0.538f)),
			};

			string path = Path.Combine(Path.GetTempPath(), $"roundtrip_mz_{Guid.NewGuid():N}.speclib");
			try
			{
				DiaNNSpecLibWriter.Write(entries, path);
				var readBack = DiaNNSpecLibReader.Read(path);

				Assert.That(readBack, Has.Count.EqualTo(1));
				var e = readBack[0];

				Assert.That((float)e.PrecursorMz, Is.EqualTo(expectedMz).Within(0.001f),
					"PrecursorMz must survive binary round-trip within float tolerance.");
				Assert.That((float)e.RetentionTime, Is.EqualTo(expectedIRT).Within(0.01f),
					"iRT must survive binary round-trip within float tolerance.");
			}
			finally
			{
				if (File.Exists(path)) File.Delete(path);
			}
		}

		/// <summary>
		/// WHAT: Fragment null placeholder slots (mz = 0.0) are discarded by the reader.
		/// WHY:  The writer always emits exactly 12 fragment slots; slots beyond the real
		///       fragment count are null. If the reader does not discard them, every entry
		///       gains phantom zero-m/z fragments that would corrupt spectral angle scoring
		///       in MetaMorpheus.
		/// </summary>
		[Test]
		public void RoundTrip_NullFragmentSlots_AreDiscarded()
		{
			// Write 3 real fragments; writer pads to 12 with null slots.
			var entries = new List<DiaNNLibraryEntry>
			{
				MakeEntry(0, "MSTYNQK", "MSTYNQK", 2, 436.2026f, 21.51f,
					DiaNNBinaryStructs.TerminusType.Internal, "PROT_A",
					(740.36f, 1.0f), (653.33f, 0.538f), (552.28f, 0.298f)),
			};

			string path = Path.Combine(Path.GetTempPath(), $"roundtrip_frags_{Guid.NewGuid():N}.speclib");
			try
			{
				DiaNNSpecLibWriter.Write(entries, path);
				var readBack = DiaNNSpecLibReader.Read(path);

				Assert.That(readBack[0].Fragments, Has.Count.EqualTo(3),
					"Reader must discard null placeholder fragment slots (mz = 0.0). " +
					"Only the 3 real fragments should survive.");
			}
			finally
			{
				if (File.Exists(path)) File.Delete(path);
			}
		}

		/// <summary>
		/// WHAT: Fragment intensities are normalized so the maximum = 1.0 on write, and
		///       the relative ratios are preserved within float precision on read-back.
		/// WHY:  DIA-NN stores pre-normalized intensities. If the writer skips normalization,
		///       all downstream spectral similarity calculations (e.g. spectral angle) produce
		///       incorrect results when the input library has unnormalized intensities.
		/// </summary>
		[Test]
		public void RoundTrip_FragmentIntensities_NormalizedAndPreserved()
		{
			// Input: max intensity is 2.0 (unnormalized), ratio between fragments is 2:1.
			var entries = new List<DiaNNLibraryEntry>
			{
				MakeEntry(0, "MSTYNQK", "MSTYNQK", 2, 436.2026f, 21.51f,
					DiaNNBinaryStructs.TerminusType.Internal, "PROT_A",
					(740.36f, 2.0f),   // should become 1.0 after normalization
                    (653.33f, 1.0f)),  // should become 0.5 after normalization
            };

			string path = Path.Combine(Path.GetTempPath(), $"roundtrip_inty_{Guid.NewGuid():N}.speclib");
			try
			{
				DiaNNSpecLibWriter.Write(entries, path);
				var readBack = DiaNNSpecLibReader.Read(path);

				var frags = readBack[0].Fragments
					.OrderByDescending(f => f.Intensity)
					.ToList();

				Assert.That(frags[0].Intensity, Is.EqualTo(1.0).Within(0.001),
					"After normalization, the highest-intensity fragment must be 1.0.");
				Assert.That(frags[1].Intensity, Is.EqualTo(0.5).Within(0.005),
					"Intensity ratio must be preserved after normalization.");
			}
			finally
			{
				if (File.Exists(path)) File.Delete(path);
			}
		}

		/// <summary>
		/// WHAT: The TerminusType on each entry survives the round-trip, specifically the
		///       two subnormal-float shared-peptide markers (SharedConflict, SharedSameType).
		/// WHY:  If the writer uses float literals for the subnormal markers instead of raw
		///       bytes, the bit pattern changes and the reader classifies them as unknown
		///       (returns null from ClassifyMarker), silently dropping all shared-peptide
		///       precursors. This was the most impactful bug found in the old code.
		/// </summary>
		[Test]
		public void RoundTrip_SharedPeptideMarkers_SurviveAsCorrectTerminusType()
		{
			var entries = new List<DiaNNLibraryEntry>
			{
                // SharedConflict — bytes 07 00 3D 00, a subnormal float
                MakeEntry(0, "GILDQSIVR", "GILDQSIVR", 2, 500.793f, 72.48f,
					DiaNNBinaryStructs.TerminusType.SharedConflict, "PROT_A",
					(830.47f, 1.0f), (717.39f, 0.98f)),

                // SharedSameType — bytes 07 00 39 00, a subnormal float
                MakeEntry(1, "MSTYNQK", "MSTYNQK", 2, 436.2026f, 21.51f,
					DiaNNBinaryStructs.TerminusType.SharedSameType, "PROT_B",
					(740.36f, 1.0f), (653.33f, 0.538f)),
			};

			string path = Path.Combine(Path.GetTempPath(), $"roundtrip_shared_{Guid.NewGuid():N}.speclib");
			try
			{
				DiaNNSpecLibWriter.Write(entries, path);
				var readBack = DiaNNSpecLibReader.Read(path);

				Assert.That(readBack, Has.Count.EqualTo(2),
					"Both shared-peptide precursors must be recovered — the subnormal " +
					"marker bytes must not be written as float literals.");

				var conflict = readBack.FirstOrDefault(e => e.StrippedSequence == "GILDQSIVR");
				Assert.That(conflict, Is.Not.Null, "GILDQSIVR (SharedConflict) must be present.");
				Assert.That(conflict!.TerminusType,
					Is.EqualTo(DiaNNBinaryStructs.TerminusType.SharedConflict));

				var same = readBack.FirstOrDefault(e => e.StrippedSequence == "MSTYNQK");
				Assert.That(same, Is.Not.Null, "MSTYNQK (SharedSameType) must be present.");
				Assert.That(same!.TerminusType,
					Is.EqualTo(DiaNNBinaryStructs.TerminusType.SharedSameType));
			}
			finally
			{
				if (File.Exists(path)) File.Delete(path);
			}
		}

		/// <summary>
		/// WHAT: The charge on each precursor survives the round-trip.
		/// WHY:  Charge is one of the two fields (with m/z) that identify a precursor for
		///       DIA window matching. A wrong charge would route a precursor to the wrong
		///       isolation window.
		/// </summary>
		[Test]
		public void RoundTrip_PrecursorCharge_IsPreserved()
		{
			var entries = new List<DiaNNLibraryEntry>
			{
				MakeEntry(0, "MSTYNQK", "MSTYNQK", 3, 291.14f, 21.51f,
					DiaNNBinaryStructs.TerminusType.Internal, "PROT_A",
					(389.21f, 1.0f)),
			};

			string path = Path.Combine(Path.GetTempPath(), $"roundtrip_charge_{Guid.NewGuid():N}.speclib");
			try
			{
				DiaNNSpecLibWriter.Write(entries, path);
				var readBack = DiaNNSpecLibReader.Read(path);

				Assert.That(readBack[0].PrecursorCharge, Is.EqualTo(3));
			}
			finally
			{
				if (File.Exists(path)) File.Delete(path);
			}
		}

		// ═══════════════════════════════════════════════════════════════════════════════
		// PART 4 — Reader: header validation
		// ═══════════════════════════════════════════════════════════════════════════════

		/// <summary>
		/// WHAT: Reading a file whose first int32 is not −10 throws FormatException.
		/// WHY:  The reader must fail fast on wrong-version files rather than silently
		///       misinterpreting arbitrary binary data as precursor records.
		/// </summary>
		[Test]
		public void Read_WrongVersionInt_ThrowsFormatException()
		{
			string path = Path.Combine(Path.GetTempPath(), $"bad_version_{Guid.NewGuid():N}.speclib");
			try
			{
				// Write version = 8 (old DIA-NN 2.0 era format, not supported by this reader)
				byte[] data = BitConverter.GetBytes(8);
				// Pad to a minimum plausible size
				var full = new byte[64];
				Buffer.BlockCopy(data, 0, full, 0, 4);
				File.WriteAllBytes(path, full);

				Assert.Throws<FormatException>(() => DiaNNSpecLibReader.Read(path),
					"A file with version ≠ −10 must throw FormatException.");
			}
			finally
			{
				if (File.Exists(path)) File.Delete(path);
			}
		}

		/// <summary>
		/// WHAT: Reading a non-existent file throws FileNotFoundException.
		/// WHY:  Basic API contract. Callers must not receive a NullReferenceException or
		///       an empty list for a missing file.
		/// </summary>
		[Test]
		public void Read_MissingFile_ThrowsFileNotFoundException()
		{
			Assert.Throws<FileNotFoundException>(
				() => DiaNNSpecLibReader.Read("/tmp/does_not_exist_ever.speclib"));
		}

		// ═══════════════════════════════════════════════════════════════════════════════
		// PART 5 — Writer: empty/null input validation
		// ═══════════════════════════════════════════════════════════════════════════════

		/// <summary>
		/// WHAT: Write with a null entries list throws ArgumentNullException.
		/// WHAT: Write with an empty entries list throws ArgumentException.
		/// WHY:  Producing a zero-byte or header-only file for null/empty input would
		///       silently create an unreadable file with no diagnostic.
		/// </summary>
		[Test]
		public void Write_NullOrEmptyEntries_ThrowsArgumentException()
		{
			string path = Path.Combine(Path.GetTempPath(), $"never_{Guid.NewGuid():N}.speclib");

			Assert.Throws<ArgumentNullException>(
				() => DiaNNSpecLibWriter.Write(null!, path),
				"Null entries must throw ArgumentNullException.");

			Assert.Throws<ArgumentException>(
				() => DiaNNSpecLibWriter.Write(new List<DiaNNLibraryEntry>(), path),
				"Empty entries list must throw ArgumentException.");

			Assert.That(File.Exists(path), Is.False,
				"No output file should be created when input validation fails.");
		}

		/// <summary>
		/// WHAT: A round-trip with multiple proteins preserves the stripped sequence of
		///       each entry, matching it to the protein from the protein table.
		/// WHY:  Protein identity (ProteinId / accession) is stored in the protein table,
		///       not inline in each precursor record. If protein table read-back is wrong,
		///       every precursor gets the wrong protein assignment.
		/// </summary>
		[Test]
		public void RoundTrip_MultipleProteins_EachEntryHasCorrectProteinId()
		{
			var entries = new List<DiaNNLibraryEntry>
			{
				MakeEntry(0, "MSTYNQK", "MSTYNQK", 2, 436.2026f, 21.51f,
					DiaNNBinaryStructs.TerminusType.NAndCTerminal, "sp|PROT_A|PROT_A_HUMAN",
					(740.36f, 1.0f)),

				MakeEntry(1, "GILDQSIVR", "GILDQSIVR", 2, 500.793f, 72.48f,
					DiaNNBinaryStructs.TerminusType.NTerminal, "sp|PROT_B|PROT_B_HUMAN",
					(830.47f, 1.0f)),
			};

			string path = Path.Combine(Path.GetTempPath(), $"roundtrip_prot_{Guid.NewGuid():N}.speclib");
			try
			{
				DiaNNSpecLibWriter.Write(entries, path);
				var readBack = DiaNNSpecLibReader.Read(path);

				Assert.That(readBack, Has.Count.EqualTo(2));

				var mstynqk = readBack.FirstOrDefault(e => e.StrippedSequence == "MSTYNQK");
				var gildq = readBack.FirstOrDefault(e => e.StrippedSequence == "GILDQSIVR");

				Assert.That(mstynqk, Is.Not.Null);
				Assert.That(gildq, Is.Not.Null);

				Assert.That(mstynqk!.ProteinId, Is.EqualTo("sp|PROT_A|PROT_A_HUMAN"),
					"MSTYNQK must be mapped to PROT_A after protein table round-trip.");
				Assert.That(gildq!.ProteinId, Is.EqualTo("sp|PROT_B|PROT_B_HUMAN"),
					"GILDQSIVR must be mapped to PROT_B after protein table round-trip.");
			}
			finally
			{
				if (File.Exists(path)) File.Delete(path);
			}
		}
	}
}