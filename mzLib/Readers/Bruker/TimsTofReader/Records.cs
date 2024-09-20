using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers.Bruker
{
    internal class Ms1Record
    {
        internal int PrecursorId { get; }
        internal int ScanStart { get; }
        internal int ScanEnd { get; }
        internal double ScanMedian { get; }

        public Ms1Record(int precursorId, int scanStart, int scanEnd, double scanMedian)
        {
            PrecursorId = precursorId;
            ScanStart = scanStart;
            ScanEnd = scanEnd;
            ScanMedian = scanMedian;
        }
    }

    //internal class Ms2Record
    //{
    //    List<long> FrameList { get; }
    //    internal int PrecursorId { get; }
    //    internal int ScanStart { get; }
    //    internal int ScanEnd { get; }
    //    internal double ScanMedian { get; }
    //    internal float IsolationMz { get; }
    //    internal float IsolationWidth { get; }
    //    internal float CollisionEnergy { get; }
    //    internal float MostAbundantPrecursorMz { get; }
    //    internal float PrecursorMonoisotopicMz { get; }
    //    internal int Charge { get; }
    //    internal float Precurso


    //    var frameList = sqliteReader.GetString(0).Split(',').Select(id => Int64.Parse(id));
    //    allFrames.UnionWith(frameList);
    //                var scanStart = sqliteReader.GetInt32(1);
    //    var scanEnd = sqliteReader.GetInt32(2);
    //    var isolationMz = sqliteReader.GetFloat(3);
    //    var isolationWidth = sqliteReader.GetFloat(4);
    //    var collisionEnergy = sqliteReader.GetFloat(5);
    //    var mostAbundantPrecursorPeak = sqliteReader.GetFloat(6);
    //    var precursorMonoisotopicMz = sqliteReader.GetFloat(7);
    //    var charge = sqliteReader.GetInt32(8);
    //    var precursorIntensity = sqliteReader.GetFloat(9);
    //    var scanMedian = sqliteReader.GetFloat(10);
    //    var precursorId = sqliteReader.GetInt32(11);
    //    runningTotal++;
    //}
}
