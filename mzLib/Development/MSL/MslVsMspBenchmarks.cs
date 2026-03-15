// Development/MSL/MslVsMspBenchmarks.cs
// Prompt 8 — Format comparison benchmarks: MSL binary vs MSP text.
//
// Run from the Development project root (NOT via dotnet test):
//   dotnet run -c Release -- --filter "*MslVsMspBenchmarks*"
//
// Requires BenchmarkDotNet in Development.csproj:
//   <PackageReference Include="BenchmarkDotNet" Version="0.14.0" />

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Text;
using BenchmarkDotNet.Attributes;
using BenchmarkDotNet.Jobs;
using MassSpectrometry;                          // ProductType, DissociationType
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;   // MslLibraryEntry, MslFragmentIon, MslFormat
using Omics.SpectrumMatch;                       // LibrarySpectrum
using Readers.SpectralLibrary;                   // MslLibrary, MslWriter, SpectralLibrary

namespace Development.MSL
{
	/// <summary>
	/// Compares load and lookup throughput between the MSL binary format and the MSP text
	/// format for the same set of synthetic precursor entries.
	///
	/// <para>
	/// <b>Baseline</b>: <see cref="MspLoad_FullIndex"/> uses <see cref="SpectralLibrary"/>,
	/// which builds a byte-offset index over the MSP file on load. All MSL benchmarks are
	/// compared against this baseline in the BenchmarkDotNet report.
	/// </para>
	///
	/// <para>
	/// <b>MSP file authorship</b>: written by iterating <see cref="LibrarySpectrum.ToString"/>
	/// via a <see cref="StreamWriter"/> — the same path MetaMorpheus uses when serialising a
	/// spectral library to disk. The MSP file is therefore a faithful ground-truth export of
	/// the same entries that are written to the .msl file.
	/// </para>
	///
	/// Attributes
	/// ----------
	/// [MemoryDiagnoser]  — reports Gen0/Gen1/Gen2 GC counts and allocated bytes per op.
	/// [SimpleJob(...)]   — targets .NET 8 with default warmup + measurement counts.
	/// </summary>
	[MemoryDiagnoser]
	[SimpleJob(RuntimeMoniker.Net80)]
	public class MslVsMspBenchmarks
	{
		// ── Parameters ───────────────────────────────────────────────────────────────

		/// <summary>
		/// Number of synthetic precursor entries written to both the MSP and MSL files.
		/// BenchmarkDotNet runs the full suite once per value.
		/// The 50 K value matches the primary comparison target in the prompt spec;
		/// 1 K and 10 K provide the low-end curve.
		/// </summary>
		[Params(1_000, 10_000, 50_000)]
		public int NPrecursors;

		// ── Private state (populated by GlobalSetup) ─────────────────────────────────

		/// <summary>
		/// Absolute path of the temporary .msl file written in <see cref="GlobalSetup"/>.
		/// Deleted in <see cref="GlobalCleanup"/>.
		/// </summary>
		private string _tempMslPath = string.Empty;

		/// <summary>
		/// Absolute path of the temporary .msp text file written in <see cref="GlobalSetup"/>.
		/// Written by calling <see cref="LibrarySpectrum.ToString"/> per entry, matching the
		/// MetaMorpheus serialisation convention.
		/// Deleted in <see cref="GlobalCleanup"/>.
		/// </summary>
		private string _tempMspPath = string.Empty;

		/// <summary>
		/// Synthetic <see cref="MslLibraryEntry"/> list generated once in
		/// <see cref="GlobalSetup"/> and used to produce both the MSL and MSP files.
		/// Also used by the lookup benchmarks to supply alternating hit/miss keys.
		/// </summary>
		private List<MslLibraryEntry> _syntheticEntries = new();

		/// <summary>
		/// Pre-converted list of <see cref="LibrarySpectrum"/> objects derived from
		/// <see cref="_syntheticEntries"/> in <see cref="GlobalSetup"/>. Used as the
		/// canonical source for writing the MSP file and for supplying lookup keys to
		/// the MSP lookup benchmark.
		/// </summary>
		private List<LibrarySpectrum> _syntheticSpectra = new();

		/// <summary>
		/// Long-lived <see cref="MslLibrary"/> opened in index-only mode during
		/// <see cref="GlobalSetup"/> and kept alive for the MSL lookup benchmark.
		/// Disposed in <see cref="GlobalCleanup"/>.
		/// </summary>
		private MslLibrary _mslLib = null!;

		/// <summary>
		/// Long-lived <see cref="SpectralLibrary"/> opened over the MSP file during
		/// <see cref="GlobalSetup"/> and kept alive for the MSP lookup benchmark.
		/// Connections closed in <see cref="GlobalCleanup"/>.
		/// </summary>
		private SpectralLibrary _mspLib = null!;

		/// <summary>
		/// Pre-built array of 1 000 lookup keys (sequence/charge pairs) used by
		/// <see cref="Msp_TryGetSpectrum_1000Lookups"/> and
		/// <see cref="Msl_TryGetLibrarySpectrum_1000Lookups"/>. Alternates between
		/// hits (entries that exist) and misses (entries that do not exist) to measure
		/// realistic mixed-workload latency. Built in <see cref="GlobalSetup"/>.
		/// </summary>
		private (string sequence, int charge)[] _lookupKeys = Array.Empty<(string, int)>();

		// ── Setup / teardown ─────────────────────────────────────────────────────────

		/// <summary>
		/// Runs once per parameter combination before any benchmark iterations begin.
		///
		/// Steps
		/// -----
		/// 1. Generate <see cref="NPrecursors"/> synthetic entries and convert to spectra.
		/// 2. Write the .msl file via <see cref="MslWriter.Write"/>.
		/// 3. Write the .msp file via <see cref="LibrarySpectrum.ToString"/> + StreamWriter.
		/// 4. Open long-lived <see cref="MslLibrary"/> (index-only) and
		///    <see cref="SpectralLibrary"/> (full index) for lookup benchmarks.
		/// 5. Build the 1 000 alternating hit/miss lookup key array.
		/// 6. Print a setup summary to stdout.
		/// </summary>
		[GlobalSetup]
		public void GlobalSetup()
		{
			var sw = Stopwatch.StartNew();

			// 1. Generate synthetic entries and convert to LibrarySpectrum for MSP writing
			_syntheticEntries = GenerateSyntheticEntries(NPrecursors);
			_syntheticSpectra = new List<LibrarySpectrum>(_syntheticEntries.Count);
			foreach (var entry in _syntheticEntries)
				_syntheticSpectra.Add(entry.ToLibrarySpectrum());

			string tempDir = Path.GetTempPath();
			string tag = Guid.NewGuid().ToString("N");

			// 2. Write .msl file
			//    Signature: MslWriter.Write(string outputPath, IReadOnlyList<MslLibraryEntry>)
			_tempMslPath = Path.Combine(tempDir, $"MslVsMsp_{NPrecursors}_{tag}.msl");
			MslWriter.Write(_tempMslPath, _syntheticEntries);

			// 3. Write .msp file using LibrarySpectrum.ToString(), matching MetaMorpheus convention
			_tempMspPath = Path.Combine(tempDir, $"MslVsMsp_{NPrecursors}_{tag}.msp");
			WriteMspFile(_syntheticSpectra, _tempMspPath);

			// 4. Open long-lived libraries for lookup benchmarks
			_mslLib = MslLibrary.LoadIndexOnly(_tempMslPath);
			_mspLib = new SpectralLibrary(new List<string> { _tempMspPath });

			// 5. Build 1 000 alternating hit/miss lookup keys
			//    Even indices → real entry (hit); odd indices → fabricated key (miss).
			_lookupKeys = new (string, int)[1_000];
			for (int i = 0; i < 1_000; i++)
			{
				if (i % 2 == 0)
				{
					// Hit: use a real entry cycling through the library
					var entry = _syntheticEntries[i % _syntheticEntries.Count];
					_lookupKeys[i] = (entry.ModifiedSequence, entry.Charge);
				}
				else
				{
					// Miss: deliberately absent sequence
					_lookupKeys[i] = ("XXXXXXXXXXX", 99);
				}
			}

			sw.Stop();

			long mslBytes = new FileInfo(_tempMslPath).Length;
			long mspBytes = new FileInfo(_tempMspPath).Length;

			Console.WriteLine();
			Console.WriteLine("MslVsMsp Benchmark Setup:");
			Console.WriteLine($"  NPrecursors:    {NPrecursors:N0}");
			Console.WriteLine($"  MSL file size:  {mslBytes / 1_048_576.0:F1} MB");
			Console.WriteLine($"  MSP file size:  {mspBytes / 1_048_576.0:F1} MB");
			Console.WriteLine($"  Size ratio:     {(double)mspBytes / mslBytes:F1}× (MSP / MSL)");
			Console.WriteLine($"  Setup time:     {sw.Elapsed.TotalSeconds:F2} s");
		}

		/// <summary>
		/// Runs once per parameter combination after all iterations complete.
		/// Disposes the long-lived library instances and deletes temp files.
		/// </summary>
		[GlobalCleanup]
		public void GlobalCleanup()
		{
			_mslLib?.Dispose();
			_mspLib?.CloseConnections();

			if (File.Exists(_tempMslPath)) File.Delete(_tempMslPath);
			if (File.Exists(_tempMspPath)) File.Delete(_tempMspPath);
		}

		// ── Load benchmarks ───────────────────────────────────────────────────────────

		/// <summary>
		/// Baseline: loads the MSP text file via <see cref="SpectralLibrary"/>, which
		/// builds a byte-offset index over the entire file before returning.
		/// All MSL benchmarks are compared against this in the BenchmarkDotNet report.
		///
		/// <see cref="BenchmarkDotNet.Attributes.BenchmarkAttribute.Baseline"/> = true causes
		/// BenchmarkDotNet to display a "Ratio" column showing how many times faster/slower
		/// each other benchmark is relative to this one.
		///
		/// Connections are closed inline so that the StreamReader handle does not leak
		/// across iterations.
		/// </summary>
		[Benchmark(Baseline = true)]
		public SpectralLibrary MspLoad_FullIndex()
		{
			var lib = new SpectralLibrary(new List<string> { _tempMspPath });
			lib.CloseConnections();
			return lib;
		}

		/// <summary>
		/// Loads the MSL file with all fragments into RAM via <see cref="MslLibrary.Load"/>.
		/// Returns the library so BenchmarkDotNet cannot elide the call; the GC finalises
		/// the handle between iterations.
		///
		/// Expected to be measurably faster than <see cref="MspLoad_FullIndex"/> for all
		/// three <see cref="NPrecursors"/> values (prompt acceptance criterion).
		/// </summary>
		[Benchmark]
		public MslLibrary MslLoad_Full()
			=> MslLibrary.Load(_tempMslPath);

		/// <summary>
		/// Loads the MSL file in index-only mode via <see cref="MslLibrary.LoadIndexOnly"/>
		/// (precursor section only; fragment data stays on disk).
		/// Expected to be the fastest load benchmark at all sizes.
		/// Returns the library so BenchmarkDotNet cannot elide the call; the GC finalises
		/// the handle between iterations.
		/// </summary>
		[Benchmark]
		public MslLibrary MslLoad_IndexOnly()
			=> MslLibrary.LoadIndexOnly(_tempMslPath);

		// ── Lookup benchmarks ─────────────────────────────────────────────────────────

		/// <summary>
		/// Performs 1 000 <see cref="SpectralLibrary.TryGetSpectrum"/> lookups against the
		/// pre-loaded MSP library, alternating between hits and misses using
		/// <see cref="_lookupKeys"/>. The bool result is consumed to prevent JIT elision.
		///
		/// Provides the MSP baseline latency for single-spectrum retrieval, against which
		/// <see cref="Msl_TryGetLibrarySpectrum_1000Lookups"/> is compared.
		/// </summary>
		[Benchmark]
		public void Msp_TryGetSpectrum_1000Lookups()
		{
			for (int i = 0; i < _lookupKeys.Length; i++)
			{
				var (seq, charge) = _lookupKeys[i];
				_ = _mspLib.TryGetSpectrum(seq, charge, out _);
			}
		}

		/// <summary>
		/// Performs 1 000 <see cref="MslLibrary.TryGetLibrarySpectrum"/> lookups against the
		/// pre-loaded MSL library (index-only mode), alternating between hits and misses using
		/// <see cref="_lookupKeys"/>. The bool result is consumed to prevent JIT elision.
		///
		/// Expected to be faster than <see cref="Msp_TryGetSpectrum_1000Lookups"/> because
		/// hits are served from the LRU cache and the dictionary lookup is O(1) in both cases,
		/// but MSL avoids the StreamReader seek and text parsing overhead on hits.
		/// </summary>
		[Benchmark]
		public void Msl_TryGetLibrarySpectrum_1000Lookups()
		{
			for (int i = 0; i < _lookupKeys.Length; i++)
			{
				var (seq, charge) = _lookupKeys[i];
				_ = _mslLib.TryGetLibrarySpectrum(seq, charge, out _);
			}
		}

		// ── MSP file writer ───────────────────────────────────────────────────────────

		/// <summary>
		/// Writes a list of <see cref="LibrarySpectrum"/> objects to an MSP text file by
		/// calling <see cref="LibrarySpectrum.ToString"/> on each entry and writing the
		/// result to a <see cref="StreamWriter"/>.
		///
		/// This matches the serialisation path used by MetaMorpheus
		/// (<c>WriteSpectrumLibrary</c> in <c>MetaMorpheusTask</c>) and by
		/// <see cref="SpectralLibrary.WriteResults"/>. Using the same path ensures the MSP
		/// file is a faithful representation that <see cref="SpectralLibrary"/> can index.
		/// </summary>
		/// <param name="spectra">Spectra to serialise. Must not be null.</param>
		/// <param name="outputPath">Destination file path. Created or overwritten.</param>
		private static void WriteMspFile(List<LibrarySpectrum> spectra, string outputPath)
		{
			using var writer = new StreamWriter(outputPath, append: false, Encoding.UTF8);
			foreach (var spectrum in spectra)
				writer.WriteLine(spectrum.ToString());
		}

		// ── Synthetic data factory ────────────────────────────────────────────────────

		/// <summary>
		/// Generates a deterministic list of <paramref name="count"/> synthetic
		/// <see cref="MslLibraryEntry"/> objects with 10 alternating b/y terminal fragment
		/// ions each — the same design as <see cref="MslBenchmarks.GenerateSyntheticEntries"/>
		/// so that results are directly comparable across the two benchmark classes.
		///
		/// Property names confirmed from MslLibraryEntry source
		/// ------------------------------------------------------
		///   ModifiedSequence, StrippedSequence, PrecursorMz (double), Charge (int),
		///   Irt (double), IsDecoy, MoleculeType, DissociationType, ProteinAccession,
		///   QValue, Fragments.
		///
		/// MslFragmentIon properties used
		/// --------------------------------
		///   ProductType, FragmentNumber, Charge, Mz, Intensity,
		///   NeutralLoss (double 0.0 = none), SecondaryProductType (null),
		///   SecondaryFragmentNumber (0).
		/// </summary>
		/// <param name="count">Number of precursor entries to generate.</param>
		/// <returns>New <see cref="List{MslLibraryEntry}"/>.</returns>
		private static List<MslLibraryEntry> GenerateSyntheticEntries(int count)
		{
			const string AminoAcids = "ACDEFGHIKLMNPQRSTVWY";
			const int AvgFragments = 10;

			var entries = new List<MslLibraryEntry>(count);

			for (int i = 0; i < count; i++)
			{
				// 7-residue sequence cycling through the amino-acid alphabet
				char[] seq = new char[7];
				for (int k = 0; k < 7; k++)
					seq[k] = AminoAcids[(i + k) % AminoAcids.Length];
				string strippedSequence = new string(seq);

				double precursorMz = 400.00 + i * 0.01;
				double irt = 10.0 + i * 0.001;

				var fragments = new List<MslFragmentIon>(AvgFragments);
				for (int f = 0; f < AvgFragments; f++)
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