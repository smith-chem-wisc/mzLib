﻿using Proteomics.PSM;

namespace Readers;

public class PsmFromTsvFile : ResultFile<PsmFromTsv>, IResultFile
{
    public override SupportedFileType FileType => SupportedFileType.psmtsv;
    public override Software Software { get; set; }

    /// <summary>
    /// Constructor used to initialize from the factory method
    /// </summary>
    public PsmFromTsvFile() : base() { }

    public PsmFromTsvFile(string filePath) : base(filePath, Software.MetaMorpheus) { }

    public override void LoadResults()
    {
        Results = SpectrumMatchTsvReader.ReadPsmTsv(FilePath, out List<string> warnings);
    }

    public override void WriteResults(string outputPath) => throw new NotImplementedException();
}