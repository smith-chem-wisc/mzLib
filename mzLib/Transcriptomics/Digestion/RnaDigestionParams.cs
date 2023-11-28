using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using MassSpectrometry;

namespace Transcriptomics
{
    public class RnaDigestionParams : IDigestionParams
    {
       
        public RnaDigestionParams(string rnase = "top-down", int maxMissedCleavages = 0, int minLength = 1, 
            int maxLength = int.MaxValue, int maxModificationIsoforms = 1024, int maxMods = 2,
            FragmentationTerminus fragmentationTerminus = FragmentationTerminus.Both) 
        {
            Rnase = RnaseDictionary.Dictionary[rnase];
            MaxMissedCleavages = maxMissedCleavages;
            MinLength = minLength;
            MaxLength = maxLength;
            MaxMods = maxModificationIsoforms;
            MaxModificationIsoforms = maxModificationIsoforms;
            FragmentationTerminus = fragmentationTerminus;
        }

        public int MaxMissedCleavages { get; set; }
        public int MinLength { get; set; }
        public int MaxLength { get; set; }
        public int MaxModificationIsoforms { get; set; }
        public int MaxMods { get; set; }
        public DigestionAgent Enzyme => Rnase;
        public Rnase Rnase { get; private set; }
        public FragmentationTerminus FragmentationTerminus { get; set; }

    }
}
