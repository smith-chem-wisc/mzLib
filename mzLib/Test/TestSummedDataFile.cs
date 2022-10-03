using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics.AminoAcidPolymer;
using Proteomics; 
using Readers;
using Proteomics.ProteolyticDigestion;
namespace Test;
[ExcludeFromCodeCoverage]
public class TestSummedDataFile
{
    private MsDataScan[] _testScan;
    private MzSpectrum _mzSpectrumA;
    private MzSpectrum _ms1;
    private MzSpectrum _ms2;
    private MsDataScan[] _scans; 
    [Test]
    public void Test
}