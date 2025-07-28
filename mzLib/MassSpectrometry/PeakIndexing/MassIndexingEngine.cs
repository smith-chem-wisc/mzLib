using MathNet.Numerics.RootFinding;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    public class MassIndexingEngine: IndexingEngine<IndexedMass>
    {
        protected override int BinsPerDalton => 1;
        public int MaxMass { get; set; } = 30000;

        public MassIndexingEngine()
        {
        }

        public static MassIndexingEngine? InitializeMassIndexingEngine(MsDataScan[] scanArray, DeconvolutionParameters deconParameters)
        {
            MassIndexingEngine newEngine = new();
            if (newEngine.IndexPeaks(scanArray, deconParameters))
                return newEngine;
            return null;
        }

        public bool IndexPeaks(MsDataScan[] scanArray, DeconvolutionParameters deconParameters, MzRange mzRange = null, double minMass = 0, int minCharge = 1)
        {
            if (scanArray.IsNullOrEmpty() || scanArray.All(p => p == null))
                return false;

            IndexedPeaks = new List<IndexedMass> [MaxMass];
            ScanInfoArray = new ScanInfo[scanArray.Length];
            for (int scanIndex = 0; scanIndex < scanArray.Length; scanIndex++)
            {
                ScanInfoArray[scanIndex] = new ScanInfo(scanArray[scanIndex].OneBasedScanNumber, scanIndex, scanArray[scanIndex].RetentionTime, scanArray[scanIndex].MsnOrder);
                var envelopes = Deconvoluter.Deconvolute(scanArray[scanIndex].MassSpectrum, deconParameters, mzRange);
                foreach (var envelope in envelopes)
                {
                    if (envelope.MonoisotopicMass < minMass || envelope.Charge < minCharge)
                        continue;
                    int roundedMass = (int)Math.Round(envelope.MonoisotopicMass * BinsPerDalton, 0);
                    IndexedPeaks[roundedMass] ??= new List<IndexedMass>();
                    IndexedPeaks[roundedMass].Add(new IndexedMass(envelope, scanArray[scanIndex].RetentionTime, scanIndex, scanArray[scanIndex].MsnOrder));
                }
            }
            if (IndexedPeaks == null || IndexedPeaks.Length == 0)
                return false;
            else
                return true;
        }
    }
}
