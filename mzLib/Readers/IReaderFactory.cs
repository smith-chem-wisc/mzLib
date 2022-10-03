using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;

namespace Readers
{
    public interface IReaderFactory
    {
        MsDataFile Reader { get; }
        MsDataFile CreateReader();
    }

    

}
