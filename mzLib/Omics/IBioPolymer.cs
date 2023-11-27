using Chemistry;
using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Omics.Modifications;

namespace Omics
{
    public interface IBioPolymer 
    {
        string Name { get; }
        string BaseSequence { get; }
        int Length { get; }
        string DatabaseFilePath { get; }
        bool IsDecoy { get; }
        bool IsContaminant { get; }
        string Organism { get; }
        string Accession { get; }
        public IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; }

      //  IEnumerable<IBioPolymerWithSetMods> Digest(DigestionParams digestionParams, List<Modification> allKnownFixedModifications,
       //     List<Modification> variableModifications, List<SilacLabel> silacLabels = null, (SilacLabel startLabel, SilacLabel endLabel)? turnoverLabels = null, bool topDownTruncationSearch = false);
    }
}
