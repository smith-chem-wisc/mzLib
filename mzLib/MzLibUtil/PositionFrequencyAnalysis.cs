using System;
using System.Collections.Generic;
using FlashLFQ;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MzLibUtil
{
    public class PositionFrequencyAnalyses
    {
        public static Dictionary<int, Dictionary<string, double>> PTMOccupancy(List<FlashLFQ.Peptide.Peptide> peptides)
        {
            
            Dictionary<int, Dictionary<string, double>> occupancy;
            foreach (var peptide in peptides.Sort(p => p.sequence))
            {
                peptide.
            }
        }
    }
}
