using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ.Interfaces
{
    public interface IFlashLfqIndexingEngine : IIndexingEngine
    {
        public SpectraFileInfo SpectraFile { get; }
        /// <summary>
        /// Paramaterless call to IndexPeaks
        /// </summary>
        /// <returns></returns>
        public bool IndexPeaks();
        public void ClearIndex();
        public void SerializeIndex();
        public void DeserializeIndex();
    }
}
