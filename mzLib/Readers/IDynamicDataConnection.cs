using System;
using System.Collections.Generic;
using System.Text;
using MassSpectrometry;

namespace Readers
{
    public interface IDynamicDataConnection
    {
        public MsDataScan GetOneBasedScanFromDynamicConnection(int oneBasedScanNumber, IFilteringParams filterParams = null);
        public void CloseDynamicConnection();
        public void InitiateDynamicConnection();
    }
}
