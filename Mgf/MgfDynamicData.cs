using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using UsefulProteomicsDatabases;

namespace IO.Mgf
{
    public class MgfDynamicData : DynamicDataConnection
    {
        private StreamReader reader;
        private Dictionary<int, long> scanToByteOffset;
        private static Regex scanNumberParser = new Regex(@"(^|\s)SCANS=(.*?)($|\D)");

        public MgfDynamicData(string filepath) : base(filepath)
        {
            InitiateDynamicConnection();
        }

        public override MsDataScan GetOneBasedScanFromDynamicConnection(int scanNumber, IFilteringParams filterParams = null)
        {
            if (reader == null)
            {
                throw new MzLibException("Cannot get scan; the dynamic data connection to " + FilePath + " has been closed!");
            }

            if (scanToByteOffset.TryGetValue(scanNumber, out long byteOffset))
            {
                // seek to the byte of the scan
                reader.BaseStream.Position = byteOffset;
                reader.DiscardBufferedData();

                return Mgf.GetNextMsDataOneBasedScanFromConnection(reader, new HashSet<int>(), filterParams, scanNumber);
            }
            else
            {
                throw new MzLibException("The specified scan number: " + scanNumber + " does not exist in " + FilePath);
            }
        }

        public override void CloseDynamicConnection()
        {
            if (reader != null)
            {
                reader.Dispose();
            }
        }

        protected override void InitiateDynamicConnection()
        {
            if (!File.Exists(FilePath))
            {
                throw new FileNotFoundException();
            }

            if (Path.GetExtension(FilePath).ToUpper() != ".MGF")
            {
                throw new InvalidDataException();
            }

            Loaders.LoadElements();
            reader = new StreamReader(FilePath);

            BuildIndex();
        }

        /// <summary>
        /// Mgf files are not indexed, so the index will need to be built.
        /// This method builds the index.
        /// </summary>
        private void BuildIndex()
        {
            scanToByteOffset = new Dictionary<int, long>();
            int oneBasedScanNumber = 0;
            long currentPositionByteOffset = 0;
            long oneBasedScanByteOffset = 0;
            bool scanHasAScanNumber = false;

            while (reader.Peek() > 0)
            {
                currentPositionByteOffset = TextFileReading.GetByteOffsetAtCurrentPosition(reader);

                string line = reader.ReadLine();

                if (line.StartsWith("BEGIN IONS", StringComparison.InvariantCultureIgnoreCase))
                {
                    oneBasedScanByteOffset = currentPositionByteOffset;
                    scanHasAScanNumber = false;
                }
                else if (line.StartsWith("SCANS=", StringComparison.InvariantCultureIgnoreCase))
                {
                    scanHasAScanNumber = true;

                    Match result = scanNumberParser.Match(line);
                    var scanString = result.Groups[2].Value;
                    oneBasedScanNumber = int.Parse(scanString);
                }
                else if (line.StartsWith("END IONS", StringComparison.InvariantCultureIgnoreCase))
                {
                    if (!scanHasAScanNumber)
                    {
                        oneBasedScanNumber++;
                    }

                    if (scanToByteOffset.ContainsKey(oneBasedScanNumber))
                    {
                        throw new MzLibException("Scan number " + oneBasedScanNumber.ToString() +
                            " appeared multiple times in " + FilePath + ", which is not allowed because we assume all scan numbers are unique.");
                    }

                    scanToByteOffset.Add(oneBasedScanNumber, oneBasedScanByteOffset);
                }
            }
        }
    }
}
