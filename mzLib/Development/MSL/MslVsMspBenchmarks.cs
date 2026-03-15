using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using BenchmarkDotNet.Attributes;
using BenchmarkDotNet.Jobs;
using MassSpectrometry;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;

namespace Development.MSL
{
	/// <summary>
	/// Format comparison benchmarks: .msl vs .msp for write throughput, full-load time,
	/// index-only-load time, and single-entry lookup latency.
	///
	/// Each benchmark is parameterised over three library sizes:
	///   NPrecursors = 1 000, 10 000, 50 000
	///
	/// The MSP baseline uses <see cref="SpectralLibrary"/> (mzLib's built-in MSP reader).
	/// MSL uses the binary format. Both operate on the same synthetic entry data generated
	/// by <see cref="MslBenchmarks.GenerateEntries"/> to ensure a fair comparison.
	///
	/// <b>File isolation note:</b> <see cref="MslLoad_IndexOnly"/> uses a separate temp file
	/// (<see cref="_mslIndexOnlyPath"/>) that is never opened by any other benchmark method
	/// in this class. This prevents the Windows file-locking error that occurs when
	/// <see cref="GlobalCleanup"/> tries to delete a file that <c>MslLoad_IndexOnly</c>'s
	/// index-only <c>FileStream</c> may still hold across BenchmarkDotNet process boundaries.
	/// </summary>
	[SimpleJob(RuntimeMoniker.Net80)]
	[MemoryDiagnoser]
	[HideColumns("Job", "RntmId", "WarmupCount", "LaunchCount", "TargetCount")]
	public class MslVsMspBenchmarks
	{
		// ── Parameters ────────────────────────────────────────────────────────────

		[Params(1_000, 10_000, 50_000)]
		public int NPrecursors { get; set; }

		[Params(10)]
		public int AvgFrag { get; set; }

		// ── State ─────────────────────────────────────────────────────────────────

		/// <summary>MSL file used by MslLoad_Full, MslLookup_BySequenceCharge, and the write benchmarks.</summary>
		private string _mslPath = null!;

		/// <summary>
		/// Dedicated MSL file used exclusively by <see cref="MslLoad_IndexOnly"/>.
		/// Kept separate from <see cref="_mslPath"/> so that GlobalCleanup's
		/// <c>File.Delete(_mslPath)</c> is never racing with an open index-only FileStream.
		/// </summary>
		private string _mslIndexOnlyPath = null!;

		private string _mspPath = null!;
		private IReadOnlyList<MslLibraryEntry> _entries = null!;
		private MslLibrary _mslLib = null!;
		private SpectralLibrary _mspLib = null!;

		// ── Setup / Cleanup ───────────────────────────────────────────────────────

		[GlobalSetup]
		public void GlobalSetup()
		{
			string tmp = Path.GetTempPath();
			string uid = Guid.NewGuid().ToString("N");

			_mslPath = Path.Combine(tmp, $"vsmsp_{NPrecursors}_{uid}.msl");
			_mslIndexOnlyPath = Path.Combine(tmp, $"vsmsp_io_{NPrecursors}_{uid}.msl");
			_mspPath = Path.Combine(tmp, $"vsmsp_{NPrecursors}_{uid}.msp");

			_entries = MslBenchmarks.GenerateEntries(NPrecursors, AvgFrag);

			// Write all three files once so load benchmarks read pre-existing data.
			// _mslIndexOnlyPath is an identical copy of _mslPath; written separately so
			// MslLoad_IndexOnly never touches _mslPath.
			MslWriter.Write(_mslPath, _entries);
			MslWriter.Write(_mslIndexOnlyPath, _entries);
			WriteMsp(_entries, _mspPath);

			// Pre-load MSL and MSP for the lookup benchmarks.
			_mslLib = MslLibrary.Load(_mslPath);
			_mspLib = new SpectralLibrary(new List<string> { _mspPath });
		}

		[GlobalCleanup]
		public void GlobalCleanup()
		{
			_mslLib?.Dispose();
			_mslLib = null!;

			_mspLib?.CloseConnections();
			_mspLib = null!;

			if (File.Exists(_mslPath)) File.Delete(_mslPath);
			if (File.Exists(_mslIndexOnlyPath)) File.Delete(_mslIndexOnlyPath);
			if (File.Exists(_mspPath)) File.Delete(_mspPath);
		}

		// ── Write: MSL vs MSP ─────────────────────────────────────────────────────

		/// <summary>
		/// Time to serialise <see cref="NPrecursors"/> entries to a new .msl binary file.
		/// </summary>
		[Benchmark]
		public void MslWrite_Full()
		{
			string path = Path.Combine(Path.GetTempPath(), $"msl_w_{Guid.NewGuid():N}.msl");
			try { MslWriter.Write(path, _entries); }
			finally { if (File.Exists(path)) File.Delete(path); }
		}

		/// <summary>
		/// Time to serialise <see cref="NPrecursors"/> entries to a new .msp text file.
		/// </summary>
		[Benchmark]
		public void MspWrite_Full()
		{
			string path = Path.Combine(Path.GetTempPath(), $"msp_w_{Guid.NewGuid():N}.msp");
			try { WriteMsp(_entries, path); }
			finally { if (File.Exists(path)) File.Delete(path); }
		}

		// ── Full load: MSL vs MSP ─────────────────────────────────────────────────

		/// <summary>
		/// Time to fully deserialise an .msl binary file into an in-memory
		/// <see cref="MslLibrary"/>. Uses <see cref="_mslPath"/>.
		/// </summary>
		[Benchmark]
		public void MslLoad_Full()
		{
			using MslLibrary lib = MslLibrary.Load(_mslPath);
		}

		/// <summary>
		/// Time to open an .msp text file and build the byte-offset index in
		/// <see cref="SpectralLibrary"/>. All work happens in the constructor;
		/// <see cref="SpectralLibrary.CloseConnections"/> releases the StreamReader
		/// at the end of each iteration.
		/// </summary>
		[Benchmark]
		public void MspLoad_Full()
		{
			var lib = new SpectralLibrary(new List<string> { _mspPath });
			lib.CloseConnections();
		}

		// ── Index-only load: MSL ──────────────────────────────────────────────────

		/// <summary>
		/// Time to open an .msl file in index-only mode (reads only the precursor index
		/// section; fragment data is left on disk), then dispose (release the FileStream).
		///
		/// Uses <see cref="_mslIndexOnlyPath"/> — a dedicated copy of the library file that
		/// is never opened by any other benchmark method. This prevents the Windows
		/// file-locking error that occurred in earlier runs when <see cref="GlobalCleanup"/>
		/// attempted to delete a path that another benchmark's in-flight FileStream still held.
		///
		/// The <c>using</c> declaration disposes the library (and releases its FileStream)
		/// before the method returns, so no handles accumulate across iterations.
		///
		/// MSP has no equivalent index-only mode (the entire text file must be scanned
		/// to build the byte-offset index), so this benchmark has no MSP counterpart.
		/// </summary>
		[Benchmark]
		public void MslLoad_IndexOnly()
		{
			using MslLibrary lib = MslLibrary.LoadIndexOnly(_mslIndexOnlyPath);
		}

		// ── Single-entry lookup: MSL vs MSP ───────────────────────────────────────

		/// <summary>
		/// Time to retrieve a single entry from the pre-loaded MSL library by sequence
		/// + charge (O(1) dictionary lookup in the in-memory index).
		/// </summary>
		[Benchmark]
		public bool MslLookup_BySequenceCharge()
		{
			return _mslLib.TryGetEntry(
				_entries[0].ModifiedSequence,
				_entries[0].Charge,
				out _);
		}

		/// <summary>
		/// Time to retrieve a single entry from the pre-loaded MSP library by sequence
		/// + charge. <see cref="SpectralLibrary.TryGetSpectrum"/> performs an O(1)
		/// dictionary lookup into its byte-offset index, then seeks to the offset and
		/// reads the spectrum record from disk (with LRU buffer on subsequent calls).
		/// </summary>
		[Benchmark]
		public bool MspLookup_BySequenceCharge()
		{
			return _mspLib.TryGetSpectrum(
				_entries[0].ModifiedSequence,
				_entries[0].Charge,
				out _);
		}

		// ── MSP write helper ──────────────────────────────────────────────────────

		/// <summary>
		/// Writes <paramref name="entries"/> to an MSP text file at <paramref name="path"/>
		/// using the standard NIST MSP format read by <see cref="SpectralLibrary"/>.
		/// </summary>
		private static void WriteMsp(IReadOnlyList<MslLibraryEntry> entries, string path)
		{
			using var sw = new StreamWriter(path, append: false,
				encoding: new UTF8Encoding(encoderShouldEmitUTF8Identifier: false),
				bufferSize: 65536);

			foreach (MslLibraryEntry e in entries)
			{
				// Normalise intensities to max = 1.0 within this precursor
				float maxInty = 0f;
				foreach (MslFragmentIon frag in e.Fragments)
					if (frag.Intensity > maxInty) maxInty = frag.Intensity;
				if (maxInty <= 0f) maxInty = 1f;

				sw.WriteLine(FormattableString.Invariant(
					$"Name: {e.ModifiedSequence}/{e.Charge}"));
				sw.WriteLine(FormattableString.Invariant(
					$"MW: {e.PrecursorMz:F6}"));
				sw.WriteLine(FormattableString.Invariant(
					$"Comment: Parent={e.PrecursorMz:F6} iRT={e.Irt:F4}"));
				sw.WriteLine(FormattableString.Invariant(
					$"Num peaks: {e.Fragments.Count}"));

				foreach (MslFragmentIon frag in e.Fragments)
				{
					float normInty = frag.Intensity / maxInty;
					sw.WriteLine(FormattableString.Invariant(
						$"{frag.Mz:F6}\t{normInty:F6}\t\"{frag.ProductType}{frag.FragmentNumber}^{frag.Charge}/0ppm\""));
				}

				sw.WriteLine(); // blank line between entries
			}
		}
	}
}