using Chemistry;
using System.Collections.Generic;
using System.Linq;

namespace IO.Thermo.Deconvolution
{
    public class IsotopicPeakGroup
    {

        #region Public Fields

        public List<ThermoMzPeak> peakList = new List<ThermoMzPeak>();
        public int charge;

        #endregion Public Fields

        #region Public Constructors

        public IsotopicPeakGroup(ThermoMzPeak peak)
        {
            peakList.Add(peak);
            this.charge = peak.Charge;
        }

        #endregion Public Constructors

        #region Public Properties

        public int Count
        {
            get
            {
                return peakList.Count;
            }
        }

        public double MostIntenseMass
        {
            get
            {
                return peakList.OrderByDescending(b => b.Y).First().X.ToMass(charge);
            }
        }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            return "C: " + charge;
        }

        #endregion Public Methods

        #region Internal Methods

        internal bool AttemptToAddNextIsotopicPeak(ThermoMzPeak peak, double tol)
        {
            if (peak.Charge == charge &&
                peak.Mz.ToMass(charge) < (peakList.Last().X.ToMass(charge) + 1) + tol && peak.Mz.ToMass(charge) > (peakList.Last().X.ToMass(charge) + 1) - tol)
            {
                peakList.Add(peak);
                return true;
            }
            return false;
        }

        #endregion Internal Methods

    }
}