using System.Text.RegularExpressions;

namespace IO.Thermo
{
    public class ThermoGlobalParams
    {
        public readonly string[] InstrumentMethods;
        public readonly string PbstrInstModel;
        public readonly string PbstrInstName;
        public readonly string PbstrInstSoftwareVersion;
        public readonly int PnControllerNumber;
        public readonly int PnControllerType;
        public readonly int PnNumInstMethods;
        public readonly string FilePath;

        public readonly ManagedThermoHelperLayer.PrecursorInfo[] CouldBePrecursor;
        public readonly int[] ScanEvent;
        public readonly int[] MsOrderByScan;

        public ThermoGlobalParams(int pnNumInstMethods, string[] instrumentMethods, string pbstrInstSoftwareVersion, string pbstrInstName, string pbstrInstModel, int pnControllerType, int pnControllerNumber, ManagedThermoHelperLayer.PrecursorInfo[] couldBePrecursor, string filePath, int[] msOrderByScan)
        {
            PnNumInstMethods = pnNumInstMethods;
            InstrumentMethods = instrumentMethods;
            PbstrInstSoftwareVersion = pbstrInstSoftwareVersion;
            PbstrInstName = pbstrInstName;
            PbstrInstModel = pbstrInstModel;
            PnControllerType = pnControllerType;
            PnControllerNumber = pnControllerNumber;
            CouldBePrecursor = couldBePrecursor;
            ScanEvent = new int[couldBePrecursor.Length];
            FilePath = filePath;
            MsOrderByScan = msOrderByScan;
        }

        public bool MonoisotopicselectionEnabled
        {
            get
            {
                foreach (var yha in InstrumentMethods)
                {
                    if (Regex.IsMatch(yha, "Monoisotopic precursor selection enabled"))
                    {
                        return true;
                    }
                }
                return false;
            }
        }
    }
}