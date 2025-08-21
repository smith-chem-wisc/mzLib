using Chemistry;
using MzLibUtil;
using MassSpectrometry; 
namespace Readers
{
    public static class MsDataFileExtensions
    {
        // <summary>
        /// Extracts an ion chromatogram from the spectra file, given a mass, charge, retention time, and mass tolerance.
        /// </summary>
        public static ExtractedIonChromatogram ExtractIonChromatogram(this MsDataFile file, double neutralMass, int charge, Tolerance massTolerance, double retentionTimeInMinutes, int msOrder = 1, double retentionTimeWindowWidthInMinutes = 5)
        {
            double theorMz = neutralMass.ToMz(charge);
            double startRt = retentionTimeInMinutes - retentionTimeWindowWidthInMinutes / 2;
            double endRt = retentionTimeInMinutes + retentionTimeWindowWidthInMinutes / 2;
            List<IIndexedPeak> xicData = new();
            IEnumerable<MsDataScan> scansInRtWindow = file.GetMsScansInTimeRange(startRt, endRt);
            foreach (MsDataScan scan in scansInRtWindow.Where(p => p.MsnOrder == msOrder))
            {
                int ind = scan.MassSpectrum.GetClosestPeakIndex(theorMz);
                double expMz = scan.MassSpectrum.XArray[ind];
                if (massTolerance.Within(expMz.ToMass(charge), neutralMass))
                {
                    xicData.Add(new IndexedMassSpectralPeak(expMz, scan.MassSpectrum.YArray[ind], scan.OneBasedScanNumber - 1, scan.RetentionTime));
                }
                else
                {
                    xicData.Add(new IndexedMassSpectralPeak(expMz, 0, scan.OneBasedScanNumber - 1, scan.RetentionTime));
                }
            }
            return new ExtractedIonChromatogram(xicData);
        }

        /// <summary>
        /// Export any MsDataFile as an MzML file to a specific file location
        /// CAUTION: some software will check the NativeID for originalScan numbers
        ///     be sure to set the new NativeID in each MsDataScan if the order has been changed
        /// </summary>
        /// <param name="file"></param>
        /// <param name="destinationPath"></param>
        /// <param name="writeIndexed"></param>
        public static void ExportAsMzML(this MsDataFile file, string destinationPath, bool writeIndexed)
        {
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(file, destinationPath, writeIndexed);
        }

        /// <summary>
        /// Creates a snip of the data file, starting at the first ms1 after the start originalScan until the end originalScan. 
        /// </summary>
        public static void ExportSnipAsMzML(this MsDataFile originalFile, int startScan, int endScan)
        {
            var filePath = originalFile.FilePath;
            if (originalFile.Scans is null)
                originalFile.LoadAllStaticData();

            var scansToKeep = new List<MsDataScan>(endScan - startScan + 1);

            bool foundFirstMs1 = false;
            int scanIndex = startScan;
            while (scanIndex < endScan + 1)
            {
                var scan = originalFile.GetOneBasedScan(scanIndex);

                // Skip until we find the first MS1 originalScan
                if (!foundFirstMs1)
                {
                    if (scan.MsnOrder == 1)
                        foundFirstMs1 = true;
                    else
                    {
                        scanIndex++;
                        continue;
                    }
                }

                scansToKeep.Add(scan);
                scanIndex++;
            }

            if (!scansToKeep.Any())
                throw new IndexOutOfRangeException(
                    $"No scans found in the range {startScan} to {endScan}. Please check the scan numbers.");

            // Replace this line
            // int scanNumberAdjustment = scansToKeep[0].OneBasedScanNumber - startScan;

            // With this line to ensure the first scan is always 1
            int scanNumberAdjustment = scansToKeep[0].OneBasedScanNumber - 1;
            var originalScanNumbers = new List<(int oneBasedScanNumber, int? oneBasedPrecursorScanNumber)>(scansToKeep.Count);
            var scanNumberMap = new Dictionary<int, int>(scansToKeep.Count * 2);
            var scanLookup = new Dictionary<int, MsDataScan>(scansToKeep.Count);
            foreach (var scan in scansToKeep)
            {
                scanLookup[scan.OneBasedScanNumber] = scan;
                originalScanNumbers.Add((scan.OneBasedScanNumber, scan.OneBasedPrecursorScanNumber));

                scanNumberMap.TryAdd(scan.OneBasedScanNumber, scan.OneBasedScanNumber - scanNumberAdjustment);
                if (scan.OneBasedPrecursorScanNumber.HasValue)
                {
                    scanNumberMap.TryAdd(scan.OneBasedPrecursorScanNumber.Value, scan.OneBasedPrecursorScanNumber.Value - scanNumberAdjustment);
                }
            }

            var scansForTheNewFile = new List<MsDataScan>(scansToKeep.Count);
            foreach (var scanNumber in originalScanNumbers.OrderBy(p => p.oneBasedScanNumber))
            {
                var originalScan = scanLookup[scanNumber.oneBasedScanNumber];
                int newScanNumber = scanNumberMap[originalScan.OneBasedScanNumber];
                int? newPrecursorScanNumber = originalScan.OneBasedPrecursorScanNumber.HasValue
                    ? scanNumberMap[originalScan.OneBasedPrecursorScanNumber.Value]
                    : null;
                string newNativeId = originalScan.NativeId.Replace($"scan={originalScan.OneBasedScanNumber}",$"scan={newScanNumber}");

                var newDataScan = new MsDataScan(
                    originalScan.MassSpectrum,
                    newScanNumber,
                    originalScan.MsnOrder,
                    originalScan.IsCentroid,
                    originalScan.Polarity,
                    originalScan.RetentionTime,
                    originalScan.ScanWindowRange,
                    originalScan.ScanFilter,
                    originalScan.MzAnalyzer,
                    originalScan.TotalIonCurrent,
                    originalScan.InjectionTime,
                    originalScan.NoiseData,
                    newNativeId,
                    originalScan.SelectedIonMZ,
                    originalScan.SelectedIonChargeStateGuess,
                    originalScan.SelectedIonIntensity,
                    originalScan.IsolationMz,
                    originalScan.IsolationWidth,
                    originalScan.DissociationType,
                    newPrecursorScanNumber,
                    originalScan.SelectedIonMonoisotopicGuessMz,
                    originalScan.HcdEnergy
                );
                scansForTheNewFile.Add(newDataScan);
            }

            string outPath = Path.Combine(Path.GetDirectoryName(filePath)!,
                filePath.GetPeriodTolerantFilenameWithoutExtension() + "_snip_" + startScan + "-" + endScan + ".mzML");

            var sourceFile = new SourceFile(
                originalFile.SourceFile.NativeIdFormat,
                originalFile.SourceFile.MassSpectrometerFileFormat,
                originalFile.SourceFile.CheckSum,
                originalFile.SourceFile.FileChecksumType,
                originalFile.SourceFile.Uri,
                originalFile.SourceFile.Id,
                originalFile.SourceFile.FileName);

            var dataFile = new GenericMsDataFile(scansForTheNewFile.ToArray(), sourceFile);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(dataFile, outPath, false);
        }
    }
}