using System.Collections.Generic;


namespace MassSpectrometry
{
    public interface IExperimentalDesign
    {
        /// <summary>
        /// A dictionary that links a file name (including the extension) to an array of sample information objects.
        /// For LFQ, the array contains one ISampleInfo object per file.
        /// For isobaric quantification (e.g., TMT), the array contains multiple ISampleInfo objects per file, corresponding to each channel.
        /// The sample info array should have a 1-to-1 mapping to the entries in an ISpectralMatch.Intensities array.
        /// </summary>
        Dictionary<string, ISampleInfo[]> FileNameSampleInfoDictionary { get; }
    }
}
