using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Proteomics
{
    public class GoTerm
    {
        public string id { get; set; }
        public string description { get; set; }
        public Aspect aspect { get; set; }
    }

    public enum Aspect
    {
        molecularFunction,
        cellularComponent,
        biologicalProcess
    }
}
