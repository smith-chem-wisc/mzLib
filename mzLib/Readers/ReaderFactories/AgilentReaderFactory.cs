using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers.ReaderFactories
{
    internal class AgilentReaderFactory : BaseReaderFactory, IReaderFactory
    {
        public MsDataFile Reader { get; }

        internal AgilentReaderFactory(string filePath) : base(filePath)
        {
            FilePath = filePath;
            Reader = CreateReader();
        }

        public MsDataFile CreateReader()
        {
            return new Agilent(FilePath);
        }
    }
}
