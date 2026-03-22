using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using BenchmarkDotNet.Attributes;
using BenchmarkDotNet.Jobs;
using MassSpectrometry;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;

namespace Development.MSL
{
	/// <summary>
	/// BenchmarkDotNet suite for core MSL library operations: write, full-load,
	/// index-only-load, m/z window query, DDA lookup, index build, and RT calibration.
	///
	/// Each benchmark is parameterised over three library sizes:
	///   NPrecursors = 1 000, 50 000, 500 000
	///
	/// All benchmarks share a common <see cref="GlobalSetup"/> that writes a synthetic
	/// .msl file to a temp path and (for query/lookup/calibration benchmarks) loads it
	/// into <see cref="_queryLib"/>. <see cref="GlobalCleanup"/> disposes and deletes
	/// the temp file.
	/// </summary>
	[SimpleJob(RuntimeMoniker.Net80)]
	[MemoryDiagnoser]
	[HideColumns("Job", "RntmId", "WarmupCount", "LaunchCount", "TargetCount")]
	public class MslBenchmarks
	{
		// ── Parameters ────────────────────────────────────────────────────────────

		[Params(1_000, 50_000, 500_000)]
		public int NPrecursors { get; set; }

		[Params(10)]
		public int AvgFrag { get; set; }

		// ── State ─────────────────────────────────────────────────────────────────

		private string _tempMslPath = null!;
		private MslLibrary _queryLib = null!;
		private IReadOnlyList<MslLibraryEntry> _entries = null!;

		// ── Setup / Cleanup ───────────────────────────────────────────────────────

		[GlobalSetup]
		public void GlobalSetup()
		{
			_tempMslPath = Path.Combine(Path.GetTempPath(),
				$"msl_bench_{NPrecursors}_{Guid.NewGuid():N}.msl");
			_entries = GenerateEntries(NPrecursors, AvgFrag);
			MslWriter.Write(_tempMslPath, _entries);
			_queryLib = MslLibrary.Load(_tempMslPath);
		}

		[GlobalCleanup]
		public void GlobalCleanup()
		{
			_queryLib?.Dispose();
			_queryLib = null!;
			if (File.Exists(_tempMslPath))
				File.Delete(_tempMslPath);
		}

		// ── Write ─────────────────────────────────────────────────────────────────

		/// <summary>
		/// Measures the time to serialise <see cref="NPrecursors"/> synthetic entries to
		/// a new .msl file. Each iteration writes to a fresh temp path to prevent OS
		/// write-caching from dominating subsequent runs.
		/// </summary>
		[Benchmark]
		public void Write_MslLibrary()
		{
			string path = Path.Combine(Path.GetTempPath(),
				$"msl_write_{Guid.NewGuid():N}.msl");
			try
			{
				MslWriter.Write(path, _entries);
			}
			finally
			{
				if (File.Exists(path))
					File.Delete(path);
			}
		}

		// ── Full load ─────────────────────────────────────────────────────────────

		/// <summary>
		/// Measures the time to fully deserialise an .msl file — header, string table,
		/// protein table, all precursor records, and all fragment blocks — into an
		/// in-memory <see cref="MslLibrary"/> ready for querying.
		/// </summary>
		[Benchmark]
		public void FullLoad_MslLibrary()
		{
			using MslLibrary lib = MslLibrary.Load(_tempMslPath);
		}

		// ── Index-only load ───────────────────────────────────────────────────────

		/// <summary>
		/// Measures the end-to-end cost of opening an .msl file in index-only mode:
		/// read the header, deserialise the precursor section into the in-memory index,
		/// and release the FileStream via Dispose.
		///
		/// The return type is <c>void</c> and each iteration uses a <c>using</c>
		/// declaration so that the <see cref="MslLibrary"/> (and its underlying
		/// FileStream) is disposed before the method returns.  This prevents multiple
		/// undisposed instances from accumulating across BenchmarkDotNet iterations,
		/// which would cause a Windows file-locking error in GlobalCleanup.
		///
		/// The measured time correctly reflects the real-world cost a caller pays when
		/// preparing an index-only library for DIA window queries: open → read precursor
		/// index → dispose.
		/// </summary>
		[Benchmark]
		public void IndexOnlyLoad_MslLibrary()
		{
			using MslLibrary lib = MslLibrary.LoadIndexOnly(_tempMslPath);
		}

		// ── m/z window queries ────────────────────────────────────────────────────

		/// <summary>
		/// Measures the latency of a single <see cref="MslLibrary.QueryMzWindow"/> call
		/// over the fully-loaded library. Uses a fixed 25 Da window centred at 800 m/z.
		/// </summary>
		[Benchmark]
		public int QueryMzWindow_SingleQuery()
		{
			ReadOnlySpan<MslPrecursorIndexEntry> results =
				_queryLib.QueryMzWindow(787.5f, 812.5f);
			return results.Length;
		}

		/// <summary>
		/// Measures the throughput of 1 000 sequential <see cref="MslLibrary.QueryMzWindow"/>
		/// calls with 25 Da windows stepping by 1 Da from 400 m/z.
		/// Represents the hot path during a DIA isolation-window sweep.
		/// </summary>
		[Benchmark]
		public long QueryMzWindow_1000Queries()
		{
			long total = 0;
			for (int i = 0; i < 1000; i++)
			{
				float centre = 400f + i;
				ReadOnlySpan<MslPrecursorIndexEntry> results =
					_queryLib.QueryMzWindow(centre - 12.5f, centre + 12.5f);
				total += results.Length;
			}
			return total;
		}

		// ── DDA lookup ────────────────────────────────────────────────────────────

		/// <summary>
		/// Measures O(1) dictionary lookup for a sequence/charge key that exists in the
		/// library. Return value prevents dead-code elimination.
		/// </summary>
		[Benchmark]
		public bool DdaLookup_ExistingEntry()
		{
			return _queryLib.TryGetEntry(
				_entries[0].FullSequence,
				_entries[0].ChargeState,
				out _);
		}

		/// <summary>
		/// Measures O(1) dictionary lookup for a sequence/charge key that does NOT exist
		/// in the library. Confirms that the negative path is equally fast.
		/// </summary>
		[Benchmark]
		public bool DdaLookup_MissingEntry()
		{
			return _queryLib.TryGetEntry("ZZZZZZNOTPRESENT", 99, out _);
		}

		// ── Index build ───────────────────────────────────────────────────────────

		/// <summary>
		/// Measures the time to construct an <see cref="MslIndex"/> from a list of
		/// already-loaded <see cref="MslLibraryEntry"/> objects (in-memory sort + index
		/// build, no I/O). The loader delegate provides direct array access so the
		/// sequence/charge dictionary is fully populated — matching the real usage pattern
		/// inside <see cref="MslLibrary.Load"/>.
		/// </summary>
		[Benchmark]
		public MslIndex BuildIndex_FromEntries()
		{
			return MslIndex.Build(_entries, i =>
				i >= 0 && i < _entries.Count ? _entries[i] : null);
		}

		// ── RT calibration ────────────────────────────────────────────────────────

		/// <summary>
		/// Measures the time to apply a linear iRT → run-RT calibration to all precursors
		/// in the loaded library, producing a new calibrated <see cref="MslLibrary"/>.
		/// This involves copying and transforming every <see cref="MslPrecursorIndexEntry"/>
		/// in the sorted index array and rebuilding the sequence/charge dictionary.
		/// </summary>
		[Benchmark]
		public MslLibrary RtCalibration_LinearTransform()
		{
			return _queryLib.WithCalibratedRetentionTimes(slope: 1.2, intercept: -5.0);
		}

		// ── Synthetic data generation ─────────────────────────────────────────────

		/// <summary>
		/// Builds a deterministic list of synthetic <see cref="MslLibraryEntry"/> objects.
		/// Precursor m/z values are spread across 400–1000 Thomson so that range queries
		/// return a realistic non-zero hit count at all three library sizes.
		/// Internal to allow reuse from <see cref="MslVsMspBenchmarks"/>.
		/// </summary>
		internal static IReadOnlyList<MslLibraryEntry> GenerateEntries(int nPrecursors, int avgFrag)
		{
			var rng = new Random(42);
			var entries = new List<MslLibraryEntry>(nPrecursors);

			string[] aa = { "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
							 "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y" };

			for (int i = 0; i < nPrecursors; i++)
			{
				// Unique sequence of 7–20 residues, deterministic from index
				int seqLen = 7 + (i % 14);
				var sb = new StringBuilder(seqLen);
				for (int j = 0; j < seqLen; j++)
					sb.Append(aa[(i * 7 + j * 13) % aa.Length]);
				string seq = sb.ToString();

				int charge = 2 + (i % 3);                                     // 2, 3, or 4
				double mz = 400.0 + (i % 1000) * 0.6 + rng.NextDouble() * 0.1;
				double irt = -30.0 + i * (120.0 / nPrecursors) + rng.NextDouble();

				var fragments = new List<MslFragmentIon>(avgFrag);
				for (int f = 0; f < avgFrag; f++)
				{
					fragments.Add(new MslFragmentIon
					{
						Mz = (float)(100.0 + f * 80.0 + rng.NextDouble() * 10),
						Intensity = (float)(rng.NextDouble() + 0.01),   // always > 0
						ProductType = (f % 2 == 0) ? ProductType.y : ProductType.b,
						FragmentNumber = f + 1,
						Charge = 1,
						NeutralLoss = 0.0
					});
				}

				entries.Add(new MslLibraryEntry
				{
					FullSequence = seq,
					BaseSequence = seq,
					PrecursorMz = mz,
					ChargeState = charge,
					RetentionTime = irt,
					IsDecoy = false,
					IsProteotypic = true,
					ProteinAccession = "BENCH_PROT",
					ProteinName = "Benchmark protein",
					GeneName = "BENCH",
					Source = MslFormat.SourceType.Predicted,
					MoleculeType = MslFormat.MoleculeType.Peptide,
					DissociationType = DissociationType.HCD,
					Nce = 28,
					MatchedFragmentIons = fragments
				});
			}

			return entries;
		}
	}
}