using Readers.QuantificationResults;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using FlashLFQ;

namespace FlashLFQ.ResultsReading
{
    public class ParsedResults
    {
        public static ParsedResults GetResultsFromFile(string peaksPath)
        {
            QuantifiedPeakFile quantResults = new QuantifiedPeakFile(peaksPath);
            //var idsFromPeaks = quantResults.MakeIdentifications();
            return null;
        }

    }
}
