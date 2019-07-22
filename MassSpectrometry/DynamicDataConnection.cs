using System;
using System.Collections.Generic;
using System.Text;

namespace MassSpectrometry
{
    public abstract class DynamicDataConnection
    {
        public readonly string FilePath;

        public DynamicDataConnection(string FilePath)
        {
            this.FilePath = FilePath;
        }

        public abstract MsDataScan GetOneBasedScanFromDynamicConnection(int oneBasedScanNumber, IFilteringParams filterParams = null);

        public abstract void CloseDynamicConnection();

        protected abstract void InitiateDynamicConnection();
    }
}
