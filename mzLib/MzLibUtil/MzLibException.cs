using System;

namespace MzLibUtil
{
    [Serializable]
    public class MzLibException(string message, Exception innerException = null) 
        : Exception(message, innerException);
}