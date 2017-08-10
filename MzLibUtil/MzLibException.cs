using System;

namespace MzLibUtil
{
    [Serializable]
    public class MzLibException : Exception
    {
        #region Public Constructors

        public MzLibException(string message) : base(message)
        {
        }

        #endregion Public Constructors
    }
}