using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry
{
    public class MassIndexingEngine: IndexingEngine<IndexedMass>
    {
        /// <summary>
        /// This will replace the default BinsPerDalton value of 100 in the base class.
        /// </summary>
        protected override int BinsPerDalton => 1;
        public int MaxMass { get; set; } = 30000;

        public MassIndexingEngine()
        {
        }

        public static MassIndexingEngine? InitializeMassIndexingEngine(MsDataScan[] scanArray, DeconvolutionParameters deconParameters)
        {
            MassIndexingEngine newEngine = new();
            if (newEngine.IndexPeaks(scanArray, deconParameters))
                return newEngine;
            return null;
        }

        /// <summary>
        /// Indexes peaks from the provided scan array using the specified deconvolution parameters.
        /// Populates the IndexedPeaks and ScanInfoArray fields.
        /// </summary>
        /// <param name="scanArray">Array of MS data scans to index.</param>
        /// <param name="deconParameters">Deconvolution parameters to use.</param>
        /// <param name="mzRange">Optional m/z range filter.</param>
        /// <param name="minMass">Minimum monoisotopic mass to consider.</param>
        /// <param name="minCharge">Minimum charge state to consider.</param>
        /// <returns>True if processing completed successfully; false if input is invalid.</returns>

        public bool IndexPeaks(MsDataScan[] scanArray, DeconvolutionParameters deconParameters, MzRange mzRange = null, double minMass = 0, int minCharge = 1)
        {
            // Validate input: return false if scan array is null, empty, or all scans are null
            if (scanArray == null || scanArray.Length == 0 || scanArray.All(p => p == null))
                return false;

            // Initialize the indexed peaks array with the maximum number of mass bins
            IndexedPeaks = new List<IndexedMass>[MaxMass];
            bool anyPeaksIndexed = false;

            // Initialize the scan info array to store metadata for each scan
            ScanInfoArray = new ScanInfo[scanArray.Length];

            // Iterate through each scan in the input array
            for (int scanIndex = 0; scanIndex < scanArray.Length; scanIndex++)
            {
                var scan = scanArray[scanIndex];
                if (scan == null)
                    continue; // Skip null scans

                // Store scan metadata
                ScanInfoArray[scanIndex] = new ScanInfo(scan.OneBasedScanNumber, scanIndex, scan.RetentionTime, scan.MsnOrder);

                // Deconvolute the scan to get isotopic envelopes
                var envelopes = Deconvoluter.Deconvolute(scan.MassSpectrum, deconParameters, mzRange);

                // Iterate through each isotopic envelope
                foreach (var envelope in envelopes)
                {
                    // Filter by minimum mass and charge
                    if (envelope.MonoisotopicMass < minMass || envelope.Charge < minCharge)
                        continue;

                    // Calculate the mass bin index for the envelope
                    int roundedMass = (int)Math.Round(envelope.MonoisotopicMass * BinsPerDalton, 0);

                    // Skip if the mass bin index is out of range
                    if (roundedMass >= IndexedPeaks.Length)
                        continue;

                    // Initialize the list for this mass bin if it doesn't exist
                    IndexedPeaks[roundedMass] ??= new List<IndexedMass>();

                    // Add the indexed mass to the appropriate mass bin
                    IndexedPeaks[roundedMass].Add(new IndexedMass(envelope, scan.RetentionTime, scanIndex, scan.MsnOrder));
                    anyPeaksIndexed = true; // Mark that at least one peak has been indexed
                }
            }

            // Throw exception if no peaks were indexed
            if (!anyPeaksIndexed)
                throw new MzLibException("No peaks in the acceptable mass or charge range.");

            // Return true to indicate successful processing
            return true;
        }
    }
}
