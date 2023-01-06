using System.Diagnostics.CodeAnalysis;
using MzLibSpectralAveraging;
using NUnit.Framework;

namespace Test.AveragingTests
{
    [ExcludeFromCodeCoverage]
    public static class GeneralMiscTests
    {
        [Test]
        public static void TestOptionToString()
        {
            SpectralAveragingOptions options1 = new SpectralAveragingOptions();
            SpectralAveragingOptions options2 = new SpectralAveragingOptions(){Percentile = 0.2, RejectionType = RejectionType.PercentileClipping};
            SpectralAveragingOptions options3 = new SpectralAveragingOptions()
                { RejectionType = RejectionType.SigmaClipping, MinSigmaValue = 1, MaxSigmaValue = 1 };
            SpectralAveragingOptions options4 = new SpectralAveragingOptions()
                { RejectionType = RejectionType.AveragedSigmaClipping, MinSigmaValue = 1, MaxSigmaValue = 1 };
            SpectralAveragingOptions options5 = new SpectralAveragingOptions()
                { RejectionType = RejectionType.WinsorizedSigmaClipping, MinSigmaValue = 1, MaxSigmaValue = 1 };
            SpectralAveragingOptions options6 = new SpectralAveragingOptions()
            {
                RejectionType = RejectionType.WinsorizedSigmaClipping, MinSigmaValue = 1, MaxSigmaValue = 1,
                SpectraWeightingType = WeightingType.NormalDistribution
            };
            SpectralAveragingOptions options7 = new SpectralAveragingOptions()
            {
                RejectionType = RejectionType.WinsorizedSigmaClipping, MinSigmaValue = 1, MaxSigmaValue = 1,
                SpectraWeightingType = WeightingType.NormalDistribution, PerformNormalization = false
            };

            string opt1String = options1.ToString();
            string opt2String = options2.ToString();
            string opt3String = options3.ToString();
            string opt4String = options4.ToString();
            string opt5String = options5.ToString();
            string opt6String = options6.ToString();
            string opt7String = options7.ToString();

            Assert.That(opt1String.Equals("NoRejection_NoWeight_Normalized_BinSize-0.01"));
            Assert.That(opt2String.Equals("PercentileClipping_NoWeight_Normalized_Percentile-0.2_BinSize-0.01"));
            Assert.That(opt3String.Equals("SigmaClipping_NoWeight_Normalized_MinSigma-1_MaxSigma-1_BinSize-0.01"));
            Assert.That(opt4String.Equals("AveragedSigmaClipping_NoWeight_Normalized_MinSigma-1_MaxSigma-1_BinSize-0.01"));
            Assert.That(opt5String.Equals("WinsorizedSigmaClipping_NoWeight_Normalized_MinSigma-1_MaxSigma-1_BinSize-0.01"));
            Assert.That(opt6String.Equals("WinsorizedSigmaClipping_NormalDistribution_Normalized_MinSigma-1_MaxSigma-1_BinSize-0.01"));
            Assert.That(opt7String.Equals("WinsorizedSigmaClipping_NormalDistribution_MinSigma-1_MaxSigma-1_BinSize-0.01"));
        }

        [Test]
        public static void TestZeroStatistics()
        {
            double[] toCalc = new double[] { 0, 0, 0, 0, 0 };

        }
    }
}
