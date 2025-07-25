using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    public class IndexedMass : IIndexedPeak
    {
        public float Intensity { get; set; }
        public float RetentionTime { get; set; }
        public int ZeroBasedScanIndex { get; set; }
        public float M { get; set; }
        public int Charge { get; set; } 
        public int MsLevel { get; set; } 
        public IsotopicEnvelope IsotopicEnvelope { get; set; }

        public IndexedMass(IsotopicEnvelope envelope, double retentionTime, int zeroBasedScanIndex, int msLevel)
        {
            IsotopicEnvelope = envelope;
            Intensity = envelope.Peaks.Max(p => (float)p.intensity);
            RetentionTime = (float)retentionTime;
            ZeroBasedScanIndex = zeroBasedScanIndex;
            M = (float)envelope.MonoisotopicMass;
            Charge = envelope.Charge;
            MsLevel = msLevel;
        }
    }
}
