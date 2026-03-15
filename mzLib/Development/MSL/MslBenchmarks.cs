// Development/Benchmarks/MslBenchmarks.cs
// Prompt 8 — BenchmarkDotNet benchmark suite for the MSL binary spectral library.
//
// ── .csproj change required before this compiles ─────────────────────────────
// Add inside an <ItemGroup> in Development/Development.csproj:
//
//   <PackageReference Include="BenchmarkDotNet" Version="0.14.0" />
//
// Run benchmarks from the Development project root (NOT via dotnet test):
//   dotnet run -c Release -- --filter "*MslBenchmarks*"
// ─────────────────────────────────────────────────────────────────────────────

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using BenchmarkDotNet.Attributes;
using BenchmarkDotNet.Jobs;
using MassSpectrometry;                          // ProductType, DissociationType
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;   // MslLibraryEntry, MslFragmentIon,
												// MslIndex, MslPrecursorIndexEntry,
												// NeutralLossCode, MoleculeType
using Readers.SpectralLibrary;                   // MslLibrary, MslWriter

namespace Development.MSL
{
	/// <summary>
	/// BenchmarkDotNet suite covering write, full-load, index-only-load, query,
	/// index-build, and RT-calibration operations on <see cref="MslLibrary"/>.
	///
	/// Attributes
	/// ----------
	/// [MemoryDiagnoser]  — reports Gen0/Gen1/Gen2 GC counts and allocated bytes per op.
	/// [SimpleJob(...)]   — targets .NET 8 with default warmup + measurement counts.
	///
	/// Design note — query benchmarks
	/// --------------------------------
	/// <see cref="MslLibrary"/> does not expose a public Index property. Query benchmarks
	/// therefore keep a long-lived <see cref="MslLibrary"/> instance alive in
	/// <see cref="_queryLib"/> (opened in index-only mode) and call its public query methods
	/// directly. The instance is disposed in <see cref="GlobalCleanup"/>.
	///
	/// Index-build benchmarks work directly with <see cref="MslIndex.Build"/>, which is
	/// the public static factory on the index type.
	///
	/// RT-calibration benchmarks call <see cref="MslLibrary.WithCalibratedRetentionTimes"/>,
	/// which returns a new <see cref="MslLibrary"/> and does not modify the original.
	/// </summary>
	[MemoryDiagnoser]
	[SimpleJob(RuntimeMoniker.Net80)]
	public class MslBenchmarks
	{
		// ── Parameters ───────────────────────────────────────────────────────────────

		/// <summary>
		/// Number of synthetic precursor entries written/read in each benchmark.
		/// BenchmarkDotNet runs the full suite once per value.
		/// Matches the three sizes in the prompt performance-target table.
		/// </summary>
		[Params(1_000, 50_000, 500_000)]
		public int NPrecursors;

		/// <summary>
		/// Average fragment ions per precursor. Fixed at 10 per the prompt spec.
		/// Declared as [Params] so the value appears in BenchmarkDotNet result tables.
		/// </summary>
		[Params(10)]
		public int AvgFragmentsPerPrecursor;

		// ── Private state (populated by GlobalSetup) ─────────────────────────────────

		/// <summary>
		/// Absolute path of the temporary .msl file written in <see cref="GlobalSetup"/>.
		/// Shared by all read/query benchmarks. Deleted in <see cref="GlobalCleanup"/>.
		/// </summary>
		private string _tempMslPath = string.Empty;

		/// <summary>
		/// Synthetic entries generated once in <see cref="GlobalSetup"/> and reused by
		/// <see cref="Write_MslLibrary"/> and <see cref="BuildIndex_FromEntries"/> so that
		/// data-generation cost is excluded from timed benchmark iterations.
		/// </summary>
		private List<MslLibraryEntry> _syntheticEntries = new();

		/// <summary>
		/// Long-lived <see cref="MslLibrary"/> opened in index-only mode during
		/// <see cref="GlobalSetup"/> and kept alive for the duration of all query benchmarks.
		/// Disposed in <see cref="GlobalCleanup"/>.
		///
		/// <see cref="MslLibrary"/> does not expose a public Index property; all query
		/// benchmarks call its public methods (<see cref="MslLibrary.QueryMzWindow"/>,
		/// <see cref="MslLibrary.TryGetEntry"/>) directly.
		/// </summary>
		private MslLibrary _queryLib = null!;

		/// <summary>
		/// Median precursor m/z of the synthetic library, pre-computed in
		/// <see cref="GlobalSetup"/>. Centres the single-query window in
		/// <see cref="QueryMzWindow_SingleQuery"/>.
		/// </summary>
		private float _medianMz;

		/// <summary>
		/// Half-width of each m/z query window: 12.5 Da each side = 25 Da total,
		/// matching the prompt specification.
		/// </summary>
		private const float QueryHalfWidth = 12.5f;

		/// <summary>
		/// 1 000 pre-computed, non-overlapping query-centre m/z values used by
		/// <see cref="QueryMzWindow_1000Queries"/>. Built in <see cref="GlobalSetup"/> by
		/// spacing 25-Da windows across the library's m/z range.
		/// </summary>
		private float[] _queryMzCentres = Array.Empty<float>();

		// ── Setup / teardown ─────────────────────────────────────────────────────────

		/// <summary>
		/// Runs once per parameter combination before any benchmark iterations begin.
		///
		/// Steps
		/// -----
		/// 1. Generate <see cref="NPrecursors"/> synthetic entries.
		/// 2. Write them to a unique temp .msl file; record wall-clock time.
		/// 3. Open the file in index-only mode and store as <see cref="_queryLib"/>.
		/// 4. Compute median m/z and 1 000 non-overlapping query window centres.
		/// 5. Print a setup summary to stdout (captured in BenchmarkDotNet logs).
		/// </summary>
		[GlobalSetup]
		public void GlobalSetup()
		{
			var sw = Stopwatch.StartNew();

			// 1. Generate entries
			_syntheticEntries = GenerateSyntheticEntries(NPrecursors, AvgFragmentsPerPrecursor);

			// 2. Write to temp file.
			//    Signature: MslWriter.Write(string outputPath, IReadOnlyList<MslLibraryEntry>)
			_tempMslPath = Path.Combine(
				Path.GetTempPath(),
				$"MslBenchmark_{NPrecursors}_{Guid.NewGuid():N}.msl");
			MslWriter.Write(_tempMslPath, _syntheticEntries);

			// 3. Open in index-only mode; keep alive for query benchmarks.
			//    MslLibrary has no public Index property — query benchmarks call its
			//    public methods (QueryMzWindow, TryGetEntry, etc.) directly.
			_queryLib = MslLibrary.LoadIndexOnly(_tempMslPath);

			// 4. Median m/z.
			//    MslLibraryEntry.PrecursorMz is double; cast to float for the array.
			var allMz = new float[_syntheticEntries.Count];
			for (int i = 0; i < _syntheticEntries.Count; i++)
				allMz[i] = (float)_syntheticEntries[i].PrecursorMz;
			Array.Sort(allMz);
			_medianMz = allMz[allMz.Length / 2];

			// 5. 1 000 non-overlapping 25-Da query windows starting below the median.
			//    Cast to float explicitly to avoid CS0266 widening to double.
			_queryMzCentres = new float[1_000];
			float start = _medianMz - (float)(1_000 * QueryHalfWidth);
			for (int i = 0; i < 1_000; i++)
				_queryMzCentres[i] = start + i * (2f * QueryHalfWidth);

			sw.Stop();

			// Heuristic size estimate: 56 B/precursor + 20 B/fragment
			long estimatedBytes =
				(long)NPrecursors * 56 +
				(long)NPrecursors * AvgFragmentsPerPrecursor * 20;
			double estimatedMb = estimatedBytes / 1_048_576.0;

			Console.WriteLine();
			Console.WriteLine("MSL Benchmark Setup:");
			Console.WriteLine($"  NPrecursors:               {NPrecursors:N0}");
			Console.WriteLine($"  AvgFragmentsPerPrecursor:  {AvgFragmentsPerPrecursor}");
			Console.WriteLine($"  Estimated file size:       {estimatedMb:F1} MB");
			Console.WriteLine($"  MSP file size (estimated): ~{estimatedMb * 10:F0} MB");
			Console.WriteLine($"  Setup time:                {sw.Elapsed.TotalSeconds:F2} s");
		}

		/// <summary>
		/// Runs once per parameter combination after all iterations complete.
		/// Disposes <see cref="_queryLib"/> (releases the index-only FileStream) and
		/// deletes the temporary .msl file.
		/// </summary>
		[GlobalCleanup]
		public void GlobalCleanup()
		{
			_queryLib?.Dispose();
			if (File.Exists(_tempMslPath))
				File.Delete(_tempMslPath);
		}

		// ── Write benchmarks ─────────────────────────────────────────────────────────

		/// <summary>
		/// Writes <see cref="_syntheticEntries"/> to a fresh unique temp file via
		/// <see cref="MslWriter.Write"/>, then immediately deletes it.
		///
		/// A new path per iteration prevents OS file-cache reuse from inflating scores.
		///
		/// Performance targets (prompt spec)
		/// ----------------------------------
		///   1 K  → &lt; 50 ms
		///  50 K  → &lt; 500 ms
		/// 500 K  → &lt; 5 s
		/// </summary>
		[Benchmark]
		public void Write_MslLibrary()
		{
			string path = Path.Combine(
				Path.GetTempPath(),
				$"MslBench_Write_{Guid.NewGuid():N}.msl");
			try
			{
				MslWriter.Write(path, _syntheticEntries);
			}
			finally
			{
				if (File.Exists(path)) File.Delete(path);
			}
		}

		// ── Read benchmarks ──────────────────────────────────────────────────────────

		/// <summary>
		/// Full-loads the pre-written temp file via <see cref="MslLibrary.Load"/>
		/// (all fragments into RAM). Returns the library so BenchmarkDotNet cannot elide
		/// the call. The GC finalises the file handle between iterations.
		///
		/// Performance targets
		/// -------------------
		///   1 K  → &lt; 10 ms
		///  50 K  → &lt; 200 ms
		/// 500 K  → &lt; 3 s
		/// </summary>
		[Benchmark]
		public MslLibrary FullLoad_MslLibrary()
			=> MslLibrary.Load(_tempMslPath);

		/// <summary>
		/// Index-only-loads the pre-written temp file via
		/// <see cref="MslLibrary.LoadIndexOnly"/> (precursor section only; fragment data
		/// stays on disk). Returns the library so BenchmarkDotNet cannot elide the call.
		/// The GC finalises the file handle between iterations.
		///
		/// Performance target: 500 K precursors → &lt; 1 s
		/// </summary>
		[Benchmark]
		public MslLibrary IndexOnlyLoad_MslLibrary()
			=> MslLibrary.LoadIndexOnly(_tempMslPath);

		// ── Query benchmarks ─────────────────────────────────────────────────────────

		/// <summary>
		/// Executes a single 25-Da m/z window query centred on <see cref="_medianMz"/>
		/// via <see cref="MslLibrary.QueryMzWindow"/>. Returns the result span so the JIT
		/// cannot dead-code-eliminate the call.
		///
		/// <see cref="MslLibrary.QueryMzWindow"/> is the zero-allocation path (delegates
		/// to <see cref="MslIndex.QueryMzRange"/> internally and returns a
		/// <see cref="ReadOnlySpan{T}"/> over the internal sorted array).
		///
		/// Performance target: median latency &lt; 1 μs (acceptance: &lt; 10 μs).
		/// </summary>
		[Benchmark]
		public ReadOnlySpan<MslPrecursorIndexEntry> QueryMzWindow_SingleQuery()
			=> _queryLib.QueryMzWindow(
				_medianMz - QueryHalfWidth,
				_medianMz + QueryHalfWidth);

		/// <summary>
		/// Executes 1 000 sequential non-overlapping 25-Da window queries via
		/// <see cref="MslLibrary.QueryMzWindow"/>. Each result span's Length is read to
		/// prevent JIT elision. Query centres were pre-computed in
		/// <see cref="GlobalSetup"/>.
		///
		/// Performance target: total elapsed &lt; 5 ms for all 1 000 queries.
		/// </summary>
		[Benchmark]
		public void QueryMzWindow_1000Queries()
		{
			for (int i = 0; i < _queryMzCentres.Length; i++)
			{
				float lo = _queryMzCentres[i] - QueryHalfWidth;
				float hi = _queryMzCentres[i] + QueryHalfWidth;
				_ = _queryLib.QueryMzWindow(lo, hi).Length;
			}
		}

		/// <summary>
		/// Measures DDA lookup latency via <see cref="MslLibrary.TryGetEntry"/> for the
		/// first entry's modified sequence and charge (always present).
		/// Returns the bool result so the JIT cannot elide the call.
		///
		/// Performance target: &lt; 100 ns (O(1) hit).
		/// </summary>
		[Benchmark]
		public bool DdaLookup_ExistingEntry()
		{
			var first = _syntheticEntries[0];
			return _queryLib.TryGetEntry(first.ModifiedSequence, first.Charge, out _);
		}

		/// <summary>
		/// Measures DDA lookup latency for a deliberately absent key via
		/// <see cref="MslLibrary.TryGetEntry"/>. Returns the bool result so the JIT
		/// cannot elide the call.
		///
		/// Performance target: &lt; 100 ns (O(1) miss).
		/// </summary>
		[Benchmark]
		public bool DdaLookup_MissingEntry()
			=> _queryLib.TryGetEntry("XXXXXXXXXXX", 99, out _);

		// ── Index-build benchmark ─────────────────────────────────────────────────────

		/// <summary>
		/// Builds an <see cref="MslIndex"/> directly from <see cref="_syntheticEntries"/>
		/// via <see cref="MslIndex.Build"/>, isolating index-construction cost from I/O.
		///
		/// Signature: Build(IReadOnlyList&lt;MslLibraryEntry&gt;, Func&lt;int, MslLibraryEntry?&gt;)
		/// The loader delegate is a closure over a local copy of the list reference so it
		/// does not close over `this`.
		///
		/// The returned index is discarded and GC-collected between iterations.
		/// </summary>
		[Benchmark]
		public MslIndex BuildIndex_FromEntries()
		{
			var entries = _syntheticEntries;
			return MslIndex.Build(
				entries,
				idx => (idx >= 0 && idx < entries.Count) ? entries[idx] : null);
		}

		// ── RT-calibration benchmark ──────────────────────────────────────────────────

		/// <summary>
		/// Applies a linear RT calibration (slope = 0.03, intercept = 40) via
		/// <see cref="MslLibrary.WithCalibratedRetentionTimes"/>.
		///
		/// This method returns a brand-new <see cref="MslLibrary"/> whose index entries
		/// have transformed Irt values; <see cref="_queryLib"/> is not modified. The new
		/// library is disposed immediately after being returned to prevent handle leaks.
		///
		/// Returning the new library prevents JIT elision; BenchmarkDotNet discards the
		/// return value automatically.
		/// </summary>
		[Benchmark]
		public MslLibrary RtCalibration_LinearTransform()
			=> _queryLib.WithCalibratedRetentionTimes(slope: 0.03, intercept: 40.0);

		// ── Synthetic data factory ────────────────────────────────────────────────────

		/// <summary>
		/// Generates a deterministic list of <paramref name="count"/> synthetic
		/// <see cref="MslLibraryEntry"/> objects with <paramref name="avgFragments"/>
		/// terminal b/y-ion pairs each.
		///
		/// Property names (confirmed from MslLibraryEntry source)
		/// -------------------------------------------------------
		///   ModifiedSequence  — string
		///   StrippedSequence  — string
		///   PrecursorMz       — double
		///   Charge            — int
		///   Irt               — double
		///   IsDecoy           — bool
		///   MoleculeType      — MslFormat.MoleculeType
		///   DissociationType  — DissociationType
		///   ProteinAccession  — string
		///   QValue            — float
		///   Fragments         — List&lt;MslFragmentIon&gt;
		///
		/// MslFragmentIon properties used
		/// --------------------------------
		///   ProductType, FragmentNumber, Charge, Mz, Intensity,
		///   NeutralLoss (double, 0.0 = none), SecondaryProductType (null),
		///   SecondaryFragmentNumber (0)
		///
		/// Design choices
		/// --------------
		/// • 7-residue sequences cycling through the 20 canonical amino acids.
		/// • Precursor m/z spaced 0.01 Da apart from 400.00 Da.
		/// • Alternating b/y ions with per-entry m/z offsets.
		/// • All entries: non-decoy, HCD, charge 2, MoleculeType.Peptide.
		/// </summary>
		/// <param name="count">Number of precursor entries to generate.</param>
		/// <param name="avgFragments">Number of fragment ions per entry.</param>
		/// <returns>New <see cref="List{MslLibraryEntry}"/>.</returns>
		private static List<MslLibraryEntry> GenerateSyntheticEntries(int count, int avgFragments)
		{
			const string AminoAcids = "ACDEFGHIKLMNPQRSTVWY";

			var entries = new List<MslLibraryEntry>(count);

			for (int i = 0; i < count; i++)
			{
				char[] seq = new char[7];
				for (int k = 0; k < 7; k++)
					seq[k] = AminoAcids[(i + k) % AminoAcids.Length];
				string strippedSequence = new string(seq);

				double precursorMz = 400.00 + i * 0.01;
				double irt = 10.0 + i * 0.001;

				var fragments = new List<MslFragmentIon>(avgFragments);
				for (int f = 0; f < avgFragments; f++)
				{
					bool isY = (f % 2) == 1;
					var productType = isY ? ProductType.y : ProductType.b;
					int fragNum = 2 + f / 2;
					float fragMz = (isY ? 200f : 150f) + fragNum * 10f + i * 0.001f;
					float intensity = Math.Max(0.1f, 1.0f - f * 0.09f);

					fragments.Add(new MslFragmentIon
					{
						ProductType = productType,
						FragmentNumber = fragNum,
						Charge = 1,
						Mz = fragMz,
						Intensity = intensity,
						NeutralLoss = 0.0,
						SecondaryProductType = null,
						SecondaryFragmentNumber = 0
					});
				}

				entries.Add(new MslLibraryEntry
				{
					ModifiedSequence = strippedSequence,
					StrippedSequence = strippedSequence,
					PrecursorMz = precursorMz,
					Charge = 2,
					Irt = irt,
					IsDecoy = false,
					MoleculeType = MslFormat.MoleculeType.Peptide,
					DissociationType = DissociationType.HCD,
					ProteinAccession = string.Empty,
					QValue = float.NaN,
					Fragments = fragments
				});
			}

			return entries;
		}
	}
}
