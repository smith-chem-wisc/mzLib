using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;

namespace Test.MslSpectralLibrary;

/// <summary>
/// <see cref="MslLibrary.LoadIndexOnly(string, int)"/> lets the caller size the entry cache of an
/// index-only library, or turn it off with 0.
/// </summary>
[TestFixture]
public sealed class TestMslIndexOnlyCacheSize
{
	private static readonly string OutputDir =
		Path.Combine(Path.GetTempPath(), "MslIndexOnlyCacheSizeTests");

	private static string Tmp(string name) => Path.Combine(OutputDir, name + ".msl");

	[OneTimeSetUp]
	public void OneTimeSetUp() => Directory.CreateDirectory(OutputDir);

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(OutputDir))
			Directory.Delete(OutputDir, recursive: true);
	}

	private static string WriteLibrary(string name, int count, int compressionLevel = 0)
	{
		var entries = new List<MslLibraryEntry>();
		for (int i = 0; i < count; i++)
		{
			entries.Add(new MslLibraryEntry
			{
				FullSequence = "PEPTIDE" + new string('K', i + 1),
				BaseSequence = "PEPTIDE" + new string('K', i + 1),
				PrecursorMz = 400.0 + 10 * i,
				ChargeState = 2,
				RetentionTime = 20.0 + i,
				Source = MslFormat.SourceType.Predicted,
				MoleculeType = MslFormat.MoleculeType.Peptide,
				DissociationType = DissociationType.HCD,
				MatchedFragmentIons = new List<MslFragmentIon>
				{
					new() { Mz = 147.11f, Intensity = 1.0f, ProductType = ProductType.y, FragmentNumber = 1, Charge = 1 },
					new() { Mz = 98.06f, Intensity = 0.5f, ProductType = ProductType.b, FragmentNumber = 2, Charge = 1 }
				}
			});
		}
		string path = Tmp(name);
		MslWriter.Write(path, entries, compressionLevel);
		return path;
	}

	[Test]
	public void OneArgumentOverload_UsesTheDefaultCache()
	{
		string path = WriteLibrary("default", 3);
		using MslLibrary lib = MslLibrary.LoadIndexOnly(path);

		Assert.That(lib.IsIndexOnly, Is.True);
		lib.GetEntry(0);
		lib.GetEntry(0);

		Assert.That(lib.IndexStatisticsForTesting.LruMisses, Is.EqualTo(1));
		Assert.That(lib.IndexStatisticsForTesting.LruHits, Is.EqualTo(1));
		Assert.That(MslLibrary.DefaultIndexOnlyCacheSize, Is.EqualTo(10_000));
	}

	[Test]
	public void ZeroCacheSize_ReadsFromDiskEveryTime_AndKeepsNothing()
	{
		string path = WriteLibrary("nocache", 3);
		using MslLibrary lib = MslLibrary.LoadIndexOnly(path, maxCachedEntries: 0);

		Assert.That(lib.IsIndexOnly, Is.True);
		MslLibraryEntry? first = lib.GetEntry(1);
		MslLibraryEntry? second = lib.GetEntry(1);

		Assert.That(first, Is.Not.Null);
		Assert.That(second, Is.Not.Null);
		Assert.That(second, Is.Not.SameAs(first), "each read builds a new entry when nothing is cached");
		Assert.That(second!.MatchedFragmentIons, Has.Count.EqualTo(2));
		Assert.That(lib.IndexStatisticsForTesting.LruHits, Is.Zero);
		Assert.That(lib.IndexStatisticsForTesting.LruMisses, Is.Zero);
	}

	[Test]
	public void SmallCache_EvictsTheEntryAddedFirst()
	{
		string path = WriteLibrary("small", 3);
		using MslLibrary lib = MslLibrary.LoadIndexOnly(path, maxCachedEntries: 2);

		lib.GetEntry(0);          // miss, cached {0}
		lib.GetEntry(1);          // miss, cached {0, 1}
		lib.GetEntry(0);          // hit
		lib.GetEntry(2);          // miss, evicts 0 (added first), cached {1, 2}
		lib.GetEntry(1);          // hit
		lib.GetEntry(0);          // miss again: 0 was evicted

		Assert.That(lib.IndexStatisticsForTesting.LruHits, Is.EqualTo(2));
		Assert.That(lib.IndexStatisticsForTesting.LruMisses, Is.EqualTo(4));
	}

	[Test]
	public void NegativeCacheSize_Throws()
	{
		string path = WriteLibrary("negative", 1);
		Assert.Throws<ArgumentOutOfRangeException>(() => MslLibrary.LoadIndexOnly(path, -1));
	}

	[Test]
	public void NullPath_Throws()
	{
		Assert.Throws<ArgumentNullException>(() => MslLibrary.LoadIndexOnly(null!, 10));
	}

	[Test]
	public void CompressedFile_LoadsInFull_WhateverTheCacheSize()
	{
		string path = WriteLibrary("compressed", 3, compressionLevel: 3);
		using MslLibrary lib = MslLibrary.LoadIndexOnly(path, maxCachedEntries: 5);

		Assert.That(lib.IsIndexOnly, Is.False);
		Assert.That(lib.IsCompressed, Is.True);
		lib.GetEntry(0);
		lib.GetEntry(0);
		Assert.That(lib.IndexStatisticsForTesting.LruHits, Is.Zero, "a fully loaded library uses no cache");
		Assert.That(lib.IndexStatisticsForTesting.LruMisses, Is.Zero);
	}

	[Test]
	public void CacheSize_SurvivesRetentionTimeCalibration()
	{
		string path = WriteLibrary("calibrated", 3);
		using MslLibrary lib = MslLibrary.LoadIndexOnly(path, maxCachedEntries: 0);
		using MslLibrary calibrated = lib.WithCalibratedRetentionTimes(1.0, 0.5);

		calibrated.GetEntry(0);
		calibrated.GetEntry(0);
		Assert.That(calibrated.IndexStatisticsForTesting.LruHits, Is.Zero);
		Assert.That(calibrated.IndexStatisticsForTesting.LruMisses, Is.Zero);
	}
}
