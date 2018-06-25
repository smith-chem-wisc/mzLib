using System.Text.RegularExpressions;

namespace IO.Thermo
{
    public class ThermoGlobalParams
    {
        public readonly string[] instrumentMethods;
        public readonly string pbstrInstModel;
        public readonly string pbstrInstName;
        public readonly string pbstrInstSoftwareVersion;
        public readonly int pnControllerNumber;
        public readonly int pnControllerType;
        public readonly int pnNumInstMethods;
        public readonly string filePath;

        public readonly ManagedThermoHelperLayer.PrecursorInfo[] couldBePrecursor;
        public readonly int[] scanEvent;
        public readonly int[] msOrderByScan;

        public ThermoGlobalParams(int pnNumInstMethods, string[] instrumentMethods, string pbstrInstSoftwareVersion, string pbstrInstName, string pbstrInstModel, int pnControllerType, int pnControllerNumber, ManagedThermoHelperLayer.PrecursorInfo[] couldBePrecursor, string filePath, int[] msOrderByScan)
        {
            this.pnNumInstMethods = pnNumInstMethods;
            this.instrumentMethods = instrumentMethods;
            this.pbstrInstSoftwareVersion = pbstrInstSoftwareVersion;
            this.pbstrInstName = pbstrInstName;
            this.pbstrInstModel = pbstrInstModel;
            this.pnControllerType = pnControllerType;
            this.pnControllerNumber = pnControllerNumber;
            this.couldBePrecursor = couldBePrecursor;
            scanEvent = new int[couldBePrecursor.Length];
            this.filePath = filePath;
            this.msOrderByScan = msOrderByScan;
        }

        public bool MonoisotopicselectionEnabled
        {
            get
            {
                foreach (var yha in instrumentMethods)
                {
                    if (Regex.IsMatch(yha, "Monoisotopic precursor selection enabled"))
                        return true;
                }
                return false;
            }
        }
    }
}