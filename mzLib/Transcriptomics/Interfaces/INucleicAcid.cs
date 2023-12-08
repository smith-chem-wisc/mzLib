using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using MassSpectrometry;
using Omics;
using Omics.Modifications;

namespace Transcriptomics
{
    public interface INucleicAcid : IHasChemicalFormula, IBioPolymer
    {
        IHasChemicalFormula FivePrimeTerminus { get; set; }

        IHasChemicalFormula ThreePrimeTerminus { get; set; }
    }
}
