using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;

namespace IO.Thermo
{
    //This file exists to deal with the configuration differences of netStandard and Net Framework, of which Thermo is only available for Net framework
    public class ThermoDataFile : MsDataFile
    {
        public ThermoGlobalParams ThermoGlobalParams { get; protected set; }

        public ThermoDataFile(int numSpectra, SourceFile sourceFile) : base(numSpectra, sourceFile)
        {

        }
        public ThermoDataFile(MsDataScan[] scans, SourceFile sourceFile) : base(scans, sourceFile)
        {

        }
    }
}
