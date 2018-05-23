using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Security.Cryptography;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace IO.Mgf
{
    public class Mgf : MsDataFile
    {
        #region Private Constructors

        private Mgf(MsDataScan[] scans, SourceFile sourceFile) : base(scans, sourceFile)
        {
        }

        #endregion Private Constructors

        #region Public Methods

        public static Mgf LoadAllStaticData(string filePath, FilteringParams filterParams = null)
        {
            if (!File.Exists(filePath))
            {
                throw new FileNotFoundException();
            }

            List<MsDataScan> scans = new List<MsDataScan>();
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
                        double precursorMass = 0;
                        int charge = 0;
                        int scanNumber = 0;
                        double rtInMinutes = 0;

                        while ((s = sr.ReadLine()) != null)
                        {
                            if (s.Equals("END IONS"))
                            {
                                readingPeaks = false;
                                MzSpectrum spectrum = new MzSpectrum(mzs.ToArray(), intensities.ToArray(), false);
                                scans.Add(new MsDataScan(spectrum, scanNumber, 2, true, charge > 0 ? Polarity.Positive : Polarity.Negative, rtInMinutes, new MzRange(mzs[0], mzs[mzs.Count]), null, MZAnalyzerType.Unknown, 0, null, null, null));
                                mzs = new List<double>();
                                intensities = new List<double>();

                                //skip the next two lines which are "" and "BEGIN IONS"
                                sr.ReadLine();
                                sr.ReadLine();
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
                                    if(sArray.Length==1)
                                    {
                                        readingPeaks = true;
                                    }
                                    else
                                    {
                                        switch(sArray[0])
                                        {
                                            case "PEPMASS":
                                                sArray = sArray[1].Split(' ');
                                                precursorMass = Convert.ToDouble(sArray[sArray.Length - 1]);
                                                break;
                                            case "CHARGE":
                                                string entry = sArray[1];
                                                charge = Convert.ToInt32(entry.Substring(0, entry.Length - 1));
                                                if(entry[entry.Length-1].Equals("-"))
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

            return new Mgf(scans.ToArray(), sourceFile);
        }

        #endregion Public Methods

    }
}
