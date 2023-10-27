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
public class TestGenericMsDataFile
{
    private MsDataScan[] _testScan;
    private MzSpectrum _mzSpectrumA;
    private MzSpectrum _ms1;
    private MzSpectrum _ms2;
    private MzSpectrum _ms12; 
    private MsDataScan[] _scans; 

    [OneTimeSetUp]
    public void OneTimeSetUp()
    {
        UsefulProteomicsDatabases.Loaders.LoadElements();

        double[] mz = { 328.73795, 329.23935, 447.73849, 448.23987, 482.23792, 482.57089, 482.90393, 500.95358, 501.28732, 501.62131, 611.99377, 612.32806, 612.66187, 722.85217, 723.35345 };
        double[] intensities = { 81007096.0, 28604418.0, 78353512.0, 39291696.0, 122781408.0, 94147520.0, 44238040.0, 71198680.0, 54184096.0, 21975364.0, 44514172.0, 43061628.0, 23599424.0, 56022696.0, 41019144.0 };

        _mzSpectrumA = new MzSpectrum(mz, intensities, false);

        var peptide = new Peptide("KQEEQMETEQQNKDEGK");

         _ms1 = CreateSpectrum(peptide.GetChemicalFormula(), 300, 2000, 1);
         _ms12 = CreateSpectrum(peptide.GetChemicalFormula(), 300, 2000, 1);
         _ms2 = CreateMS2spectrum(peptide.Fragment(FragmentTypes.b | FragmentTypes.y, true), 100, 1500);

        _scans = new MsDataScan[2];
        _scans[0] = new MsDataScan(_ms1, 1, 1, false, Polarity.Positive, 1.0, new MzRange(300, 2000), "first spectrum", MZAnalyzerType.Unknown, _ms1.SumOfAllY, 1, null, "scan=1");
        _scans[1] = new MsDataScan(_ms2, 2, 2, false, Polarity.Positive, 2.0, new MzRange(100, 1500), "second spectrum", MZAnalyzerType.Unknown, _ms2.SumOfAllY, 1, null, "scan=2", 693.9892, 3, .3872, 693.99, 1, DissociationType.Unknown, 1, 693.6550);
    }
    [Test]
    public void TestGenericDataFile()
    {
        // Tests to see if any errors are thrown in the constructors
        var sf = new SourceFile("no nativeID format", "mgf format", null, null, null);
        string dummyPath = String.Empty;
        GenericMsDataFile gFile2 = new GenericMsDataFile(_scans, sf);
        GenericMsDataFile gFile3 = new GenericMsDataFile(_scans.Length, sf);
        GenericMsDataFile gFile4 = new GenericMsDataFile(dummyPath); 
    }

    [Test]
    public void TestAbstractOverrides()
    {
        GenericMsDataFile gFile = new GenericMsDataFile("");
        Assert.Throws<NotImplementedException>(() =>
        {
            gFile.LoadAllStaticData();
        });
        Assert.Throws<NotImplementedException>(() =>
        {
            gFile.GetSourceFile();
        }); 
        Assert.Throws<NotImplementedException>(() =>
        {
            gFile.GetOneBasedScanFromDynamicConnection(1);
        }); 
        Assert.Throws<NotImplementedException>(() =>
        {
            gFile.CloseDynamicConnection();
        }); 
        Assert.Throws<NotImplementedException>(() =>
        {
            gFile.InitiateDynamicConnection();
        }); 
    }

    private MzSpectrum CreateSpectrum(ChemicalFormula f, double lowerBound, double upperBound, int minCharge)
    {
        IsotopicDistribution isodist = IsotopicDistribution.GetDistribution(f, 0.1, 0.001);
        MzSpectrum notActuallyMzS = new MzSpectrum(isodist.Masses.ToArray(), isodist.Intensities.ToArray(), false);

        notActuallyMzS.ReplaceXbyApplyingFunction(b => b.Mz.ToMz(1));

        List<double> allMasses = new List<double>();
        List<double> allIntensitiess = new List<double>();

        while (notActuallyMzS.FirstX > lowerBound)
        {
            for (int i = 0; i < notActuallyMzS.Size; i++)
            {
                if (notActuallyMzS.XArray[i] > lowerBound && notActuallyMzS.XArray[i] < upperBound)
                {
                    allMasses.Add(notActuallyMzS.XArray[i]);
                    allIntensitiess.Add(notActuallyMzS.YArray[i]);
                }
            }
            minCharge += 1;
            notActuallyMzS = new MzSpectrum(isodist.Masses.ToArray(), isodist.Intensities.ToArray(), false);
            notActuallyMzS.ReplaceXbyApplyingFunction(s => s.Mz.ToMz(minCharge));
        }
             
        var allMassesArray = allMasses.ToArray();
        var allIntensitiessArray = allIntensitiess.ToArray();

        Array.Sort(allMassesArray, allIntensitiessArray);

        return new MzSpectrum(allMassesArray, allIntensitiessArray, false);
    }
    private MzSpectrum CreateMS2spectrum(IEnumerable<Fragment> fragments, int v1, int v2)
    {
        List<double> allMasses = new List<double>();
        List<double> allIntensities = new List<double>();
        foreach (ChemicalFormulaFragment f in fragments)
        {
            var spec = CreateSpectrum(f.ThisChemicalFormula, v1, v2, 2);
            for (int i = 0; i < spec.Size; i++)
            {
                allMasses.Add(spec.XArray[i]);
                allIntensities.Add(spec.YArray[i]);
            }
        }
        var allMassesArray = allMasses.ToArray();
        var allIntensitiessArray = allIntensities.ToArray();

        Array.Sort(allMassesArray, allIntensitiessArray);
        return new MzSpectrum(allMassesArray, allIntensitiessArray, false);
    }
}