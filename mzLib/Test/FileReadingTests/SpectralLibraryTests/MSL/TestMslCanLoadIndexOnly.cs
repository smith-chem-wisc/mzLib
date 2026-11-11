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
/// <see cref="MslLibrary.CanLoadIndexOnly"/> tells a caller, before loading, whether
/// <see cref="MslLibrary.LoadIndexOnly(string)"/> will open a file index-only or fall back to a
/// full load because the file is compressed.
/// </summary>
[TestFixture]
public sealed class TestMslCanLoadIndexOnly
{
	private static readonly string OutputDir =
		Path.Combine(Path.GetTempPath(), "MslCanLoadIndexOnlyTests");

	private static string Tmp(string name) => Path.Combine(OutputDir, name);

	[OneTimeSetUp]
	public void OneTimeSetUp() => Directory.CreateDirectory(OutputDir);

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(OutputDir))
			Directory.Delete(OutputDir, recursive: true);
	}

	private static string WriteLibrary(string name, int compressionLevel)
	{
		var entries = new List<MslLibraryEntry>
		{
			new()
			{
				FullSequence = "PEPTIDEK",
				BaseSequence = "PEPTIDEK",
				PrecursorMz = 465.23,
				ChargeState = 2,
				RetentionTime = 25.0,
				Source = MslFormat.SourceType.Predicted,
				MoleculeType = MslFormat.MoleculeType.Peptide,
				DissociationType = DissociationType.HCD,
				MatchedFragmentIons = new List<MslFragmentIon>
				{
					new() { Mz = 147.11f, Intensity = 1.0f, ProductType = ProductType.y, FragmentNumber = 1, Charge = 1 },
					new() { Mz = 98.06f, Intensity = 0.5f, ProductType = ProductType.b, FragmentNumber = 2, Charge = 1 }
				}
			}
		};
		string path = Tmp(name + ".msl");
		MslWriter.Write(path, entries, compressionLevel);
		return path;
	}

	[Test]
	public void UncompressedFile_CanLoadIndexOnly_AndDoes()
	{
		string path = WriteLibrary("plain", compressionLevel: 0);

		Assert.That(MslLibrary.CanLoadIndexOnly(path), Is.True);
		using MslLibrary lib = MslLibrary.LoadIndexOnly(path);
		Assert.That(lib.IsIndexOnly, Is.True);
	}

	[Test]
	public void CompressedFile_CannotLoadIndexOnly_AndLoadIndexOnlyFallsBack()
	{
		string path = WriteLibrary("zstd", compressionLevel: 3);

		Assert.That(MslLibrary.CanLoadIndexOnly(path), Is.False);
		using MslLibrary lib = MslLibrary.LoadIndexOnly(path);
		Assert.That(lib.IsIndexOnly, Is.False);
		Assert.That(lib.IsCompressed, Is.True);
	}

	[Test]
	public void NullPath_Throws() =>
		Assert.Throws<ArgumentNullException>(() => MslLibrary.CanLoadIndexOnly(null!));

	[Test]
	public void MissingFile_Throws() =>
		Assert.Throws<FileNotFoundException>(() => MslLibrary.CanLoadIndexOnly(Tmp("missing.msl")));

	[Test]
	public void NotAnMslFile_Throws()
	{
		string path = Tmp("not_msl.msl");
		File.WriteAllText(path, "Name: PEPTIDEK/2\nNum peaks: 0\n");
		Assert.Throws<FormatException>(() => MslLibrary.CanLoadIndexOnly(path));
	}
}
