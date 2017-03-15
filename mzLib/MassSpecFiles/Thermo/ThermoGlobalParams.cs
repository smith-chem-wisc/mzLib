using System.Text.RegularExpressions;

namespace IO.Thermo
{
    public class ThermoGlobalParams
    {

        #region Public Fields

        public readonly string[] instrumentMethods;
        public readonly string pbstrInstModel;
        public readonly string pbstrInstName;
        public readonly string pbstrInstSoftwareVersion;
        public readonly int pnControllerNumber;
        public readonly int pnControllerType;
        public readonly int pnNumInstMethods;

        public readonly ClassLibrary1.PrecursorInfo[] couldBePrecursor;

        #endregion Public Fields

        #region Public Constructors

        public ThermoGlobalParams(int pnNumInstMethods, string[] instrumentMethods, string pbstrInstSoftwareVersion, string pbstrInstName, string pbstrInstModel, int pnControllerType, int pnControllerNumber, ClassLibrary1.PrecursorInfo[] couldBePrecursor)
        {
            this.pnNumInstMethods = pnNumInstMethods;
            this.instrumentMethods = instrumentMethods;
            this.pbstrInstSoftwareVersion = pbstrInstSoftwareVersion;
            this.pbstrInstName = pbstrInstName;
            this.pbstrInstModel = pbstrInstModel;
            this.pnControllerType = pnControllerType;
            this.pnControllerNumber = pnControllerNumber;
            this.couldBePrecursor = couldBePrecursor;
        }

        #endregion Public Constructors

        #region Public Properties

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

        #endregion Public Properties

    }
}