using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Transcriptomics
{
    public interface INucleicAcid
    {
        /// <summary>
        /// The amino acid sequence
        /// </summary>
        string Sequence { get; }

        /// <summary>
        /// The length of the amino acid sequence
        /// </summary>
        int Length { get; }
    }
}
