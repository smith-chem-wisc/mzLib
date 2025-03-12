using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MzLibUtil.SparseMatrix
{
    public interface ITimsTofMatrix
    {
        List<IIndexedTimsTofPeaks>[] IndexedPeaks { get; }
        int BinsPerDalton { get; }
        public void InitializeIndexedPeaksArray(MzRange ms1ScanRange);
        public IIndexedTimsTofPeaks GetIndexedPeaks(double theoreticalMass, int zeroBasedScanIndex, Tolerance tolerance, int chargeState);
    }

    public interface IIndexedTimsTofPeaks
    {
        double RetentionTime { get; }
        double Mz { get; }
        List<int> PushIndex { get; }
        List<double> Intensity { get; }
        int FrameNumber { get; }
    }
}
