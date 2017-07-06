using System;

namespace IO.MzML
{
    [Serializable]
    internal class MzmlReaderException : Exception
    {

        #region Public Constructors

        public MzmlReaderException(string message) : base(message)
        {
        }

        #endregion Public Constructors

    }
}