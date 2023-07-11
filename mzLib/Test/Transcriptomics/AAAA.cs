using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using NUnit.Framework;
using Proteomics.Fragmentation;
using Transcriptomics;

namespace Test.Transcriptomics
{
    internal class AAAA
    {

        [Test]
        public static void TryMatchingIons()
        {
            string scansPath = @"D:\Projects\RNA\TestData\SixMerMs2Scans.toml";
            RNA rna = new RNA("GUACUG");
            var scans = Toml.ReadFile<MsDataScan[]>(scansPath);
            SearchParameters parameters = new();





        }

        public List<MatchedFragmentIon> MatchFragmentIons(NucleicAcid nucleicAcid, MsDataScan scan, SearchParameters parameters)
        {




            return null;
        }
    }

    public class SearchParameters
    {
        public Tolerance ProductMassTolerance { get; set; }

        public SearchParameters(double productTolerance = 10)
        {
            ProductMassTolerance = new PpmTolerance(productTolerance);
        }

        
    }
}
