using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Transcriptomics;
using UsefulProteomicsDatabases;

namespace Readers
{
    public static class NucleicAcidDbLoader
    {
        public static IEnumerable<RNA> LoadRNAFasta(string dbLocation, bool generateTargets, DecoyType decoyType)
        {

            using var stream = new FileStream(dbLocation, FileMode.Open, FileAccess.Read, FileShare.Read)
            {

            };



            return null;
        }
    }
}
