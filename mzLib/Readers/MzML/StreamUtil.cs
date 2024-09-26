using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers
{
    public static class StreamUtil
    {
        public static void ReadExactly(Stream input, byte[] buffer, int bytesToRead)
        {
            int index = 0;
            while (index < bytesToRead)
            {
                int read = input.Read(buffer, index, bytesToRead - index);
                if (read == 0)
                {
                    throw new EndOfStreamException
                    (String.Format("End of stream reached with {0} byte{1} left to read.",
                        bytesToRead - index,
                        bytesToRead - index == 1 ? "s" : ""));
                }
                index += read;
            }
        }

    }
}
