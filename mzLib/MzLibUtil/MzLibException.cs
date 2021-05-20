using System;

namespace MzLibUtil
{
    [Serializable]
    public class MzLibException : Exception
    {
        public MzLibException(string message)
            : base(message)
        {
        }
    }
}