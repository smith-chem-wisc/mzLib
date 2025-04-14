using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers.timsTOF
{
    /// <summary>
    /// This is similar to an mz spectrum, but much more lightweight
    /// It stores intensities as ints and tof indices instead of mz values
    /// </summary>
    public class TimsSpectrum
    {
        public int[] TofIndices { get; init; }
        public int[] Intensities { get; init; }
        public int TimsIndex { get; init; }

        public TimsSpectrum(int[] tofIndices, int[] intensities, int timsIndex)
        {
            TofIndices = tofIndices;
            Intensities = intensities;
            TimsIndex = timsIndex;
        }
    }
}
