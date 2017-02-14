using System;

namespace UsefulProteomicsDatabases
{
    [Serializable]
    public class PtmListLoaderException : Exception
    {

        #region Public Constructors

        public PtmListLoaderException(string message) : base(message)
        {
        }

        #endregion Public Constructors

    }
}