using System;

namespace IO.Thermo
{
    [Serializable]
    public class ThermoReadException : Exception
    {

        #region Public Constructors

        public ThermoReadException(string message) : base(message)
        {
        }

        #endregion Public Constructors

    }
}