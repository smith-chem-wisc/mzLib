using Readers.QuantificationResults;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using FlashLFQ;
using Readers;

namespace FlashLFQ.ResultsReading
{
    public class ParsedResults
    {
        public static ParsedResults GetResultsFromFile(string peaksPath, string psmResultPath)
        {
            QuantifiedPeakFile quantResults = new QuantifiedPeakFile(peaksPath);
            PsmFromTsvFile psmResults = new PsmFromTsvFile(psmResultPath);
            var idsFromPeaks = psmResults.MakeIdentifications();

            return null;
        }

    }
}
