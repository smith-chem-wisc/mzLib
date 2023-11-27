using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Transcriptomics.Modifications
{
    public class ModomicsDbMod
    {
        public int ModomicsId { get; set; }
        public string Abbreviation { get; set; }
        public string ChemicalFormula { get; set; }
        public string Lc_Elution_Time { get; set; }
        public string AverageMass { get; set; }
        public string MonoisotopicMass { get; set; }
        public string ProtMass { get; set; }
        public string ModomicsCode { get; set; }
        public string Name { get; set; }
        public string ProductIons { get; set; }
        public string ReferenceMoiety { get; set; }
        public string ShortName { get; set; }
        public string Smile { get; set; }
    }
}
