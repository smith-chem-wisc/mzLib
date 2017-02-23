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

        public GoTerm()
        {

        }

        public GoTerm(string _id, string _description, Aspect _aspect)
        {
            this.id = _id;
            this.description = _description;
            this.aspect = _aspect;
        }
    }

    public enum Aspect
    {
        molecularFunction,
        cellularComponent,
        biologicalProcess
    }
}
