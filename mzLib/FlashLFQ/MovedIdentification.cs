using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Transactions;

namespace FlashLFQ
{
    internal class MovedIdentification : Identification
    {
        public readonly double TimeShift;
        public MovedIdentification(SpectraFileInfo fileInfo, string BaseSequence, string ModifiedSequence,
            double monoisotopicMass,
            double ms2RetentionTimeInMinutes, int chargeState, List<ProteinGroup> proteinGroups, double timeShift) : base(fileInfo, BaseSequence, ModifiedSequence,
            monoisotopicMass,ms2RetentionTimeInMinutes + timeShift, chargeState, proteinGroups ) 
        {
            
            this.TimeShift = timeShift;
        }
    }
}
