using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    public class IndexedMass : IIndexedPeak
    {
        public double Intensity { get; set; }
        public double RetentionTime { get; set; }
        public int ZeroBasedScanIndex { get; set; }
        public double M { get; set; }
        public int Charge { get; set; } 
        public int MsLevel { get; set; } 
        public IsotopicEnvelope IsotopicEnvelope { get; set; }

        public IndexedMass(IsotopicEnvelope envelope, double retentionTime, int zeroBasedScanIndex, int msLevel)
        {
            Intensity = envelope.Peaks.Max(P => P.intensity);
            RetentionTime = retentionTime;
            ZeroBasedScanIndex = zeroBasedScanIndex;
            M = envelope.MonoisotopicMass;
            Charge = envelope.Charge;
            MsLevel = msLevel;
        }
    }
}
