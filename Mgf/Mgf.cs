using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace IO.Mgf
{
    public class Mgf : MsDataFile
    {
        private MsDataScan[] indexedScans;

        private Mgf(MsDataScan[] scans, SourceFile sourceFile) : base(scans, sourceFile)
        {
            indexedScans = new MsDataScan[scans[scans.Length - 1].OneBasedScanNumber];
            foreach (MsDataScan scan in scans)
            {
                indexedScans[scan.OneBasedScanNumber - 1] = scan;
            }
        }

        public static Mgf LoadAllStaticData(string filePath, FilteringParams filterParams = null)
        {
            if (!File.Exists(filePath))
            {
                throw new FileNotFoundException();
            }

            List<MsDataScan> scans = new List<MsDataScan>();
            HashSet<int> checkForDuplicateScans = new HashSet<int>();

            using (FileStream fs = new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                using (BufferedStream bs = new BufferedStream(fs))
                {
                    using (StreamReader sr = new StreamReader(bs))
                    {
                        string s;
                        while ((s = sr.ReadLine()) != null && !s.Equals("BEGIN IONS"))
                        {
                            //do nothing with the first few scans
                        }
                        bool readingPeaks = false;
                        List<double> mzs = new List<double>();
                        List<double> intensities = new List<double>();
                        double precursorMz = 0;
                        int charge = 2;
                        int scanNumber = 0;
                        int oldScanNumber = 0;
                        double rtInMinutes = 0;

                        while ((s = sr.ReadLine()) != null)
                        {
                            if (s.Equals("END IONS"))
                            {
                                if (!checkForDuplicateScans.Add(scanNumber)) //returns false if the scan already exists
                                {
                                    throw new MzLibException("Scan number " + scanNumber.ToString() + " appeared multiple times in " + filePath + ", which is not allowed because we assume all scan numbers are unique.");
                                }

                                readingPeaks = false;
                                MzSpectrum spectrum = new MzSpectrum(mzs.ToArray(), intensities.ToArray(), false);
                                scans.Add(new MsDataScan(spectrum, scanNumber, 2, true, charge > 0 ? Polarity.Positive : Polarity.Negative, rtInMinutes, new MzRange(mzs[0], mzs[mzs.Count - 1]), null, MZAnalyzerType.Unknown, intensities.Sum(), 0, null, null, precursorMz, charge, null, precursorMz, null, DissociationType.Unknown, null, precursorMz));
                                mzs = new List<double>();
                                intensities = new List<double>();
                                oldScanNumber = scanNumber;
                                charge = 2; //default when unknown

                                //skip the next two lines which are "" and "BEGIN IONS"
                                while ((s = sr.ReadLine()) != null && !s.Equals("BEGIN IONS"))
                                {
                                    //do nothing
                                }
                            }
                            else
                            {
                                if (readingPeaks)
                                {
                                    string[] sArray = s.Split(' ');
                                    mzs.Add(Convert.ToDouble(sArray[0]));
                                    intensities.Add(Convert.ToDouble(sArray[1]));
                                }
                                else
                                {
                                    string[] sArray = s.Split('=');
                                    if (sArray.Length == 1)
                                    {
                                        readingPeaks = true;
                                        sArray = s.Split(' ');
                                        mzs.Add(Convert.ToDouble(sArray[0]));
                                        intensities.Add(Convert.ToDouble(sArray[1]));

                                        if (oldScanNumber == scanNumber) //if there's no recorded scan number, simply index them.
                                        {
                                            scanNumber++;
                                        }
                                    }
                                    else
                                    {
                                        switch (sArray[0])
                                        {
                                            case "PEPMASS":
                                                sArray = sArray[1].Split(' ');
                                                precursorMz = Convert.ToDouble(sArray[0]);
                                                break;

                                            case "CHARGE":
                                                string entry = sArray[1];
                                                charge = Convert.ToInt32(entry.Substring(0, entry.Length - 1));
                                                if (entry[entry.Length - 1].Equals("-"))
                                                {
                                                    charge *= -1;
                                                }
                                                break;

                                            case "SCANS":
                                                scanNumber = Convert.ToInt32(sArray[1]);
                                                break;

                                            case "RTINSECONDS":
                                                rtInMinutes = Convert.ToDouble(sArray[sArray.Length - 1]) / 60.0;
                                                break;

                                            default:
                                                break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            SourceFile sourceFile = new SourceFile("no nativeID format", "mgf format", null, null, null);

            return new Mgf(scans.OrderBy(x=>x.OneBasedScanNumber).ToArray(), sourceFile);
        }

        public override MsDataScan GetOneBasedScan(int scanNumber)
        {
            return indexedScans[scanNumber - 1];
        }
    }
}
