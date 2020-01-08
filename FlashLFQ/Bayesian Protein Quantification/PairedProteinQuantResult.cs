using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ
{
    public class PairedProteinQuantResult : ProteinQuantificationEngineResult
    {
        public readonly List<(Peptide peptide, List<double> foldChanges)> PeptideFoldChangeMeasurements;

        public PairedProteinQuantResult(ProteinGroup protein, string condition1, string condition2, double[] mus, double[] sds, double[] nus, 
            List<(Peptide peptide, List<double> foldChanges)> fcs) 
            : base(protein, condition1, condition2, mus, sds, nus)
        {
            this.PeptideFoldChangeMeasurements = fcs;
            //this.ConditionIntensityPointEstimate = Math.Pow(2, FoldChangePointEstimate) * referenceIntensity;
        }
    }
}
