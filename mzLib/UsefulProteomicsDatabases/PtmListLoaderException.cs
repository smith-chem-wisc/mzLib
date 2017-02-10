using System;
using System.Runtime.Serialization;

namespace UsefulProteomicsDatabases
{
    [Serializable]
    public class PtmListLoaderException : Exception
    {
        public PtmListLoaderException()
        {
        }

        public PtmListLoaderException(string message) : base(message)
        {
        }

        public PtmListLoaderException(string message, Exception innerException) : base(message, innerException)
        {
        }

        protected PtmListLoaderException(SerializationInfo info, StreamingContext context) : base(info, context)
        {
        }
    }
}