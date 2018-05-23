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

        private Mzml(MsDataScan[] scans, SourceFile sourceFile) : base(scans, sourceFile)
        {
        }

        #endregion Private Constructors

        #region Public Methods

        public static Mzml LoadAllStaticData(string filePath, FilteringParams filterParams = null)
        {
            if (!File.Exists(filePath))
            {
                throw new FileNotFoundException();
            }

            Generated.mzMLType _mzMLConnection;

            try
            {
                using (FileStream fs = new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.Read))
                {
                    var _indexedmzMLConnection = (Generated.indexedmzML)MzmlMethods.indexedSerializer.Deserialize(fs);
                    _mzMLConnection = _indexedmzMLConnection.mzML;
                }
            }
            catch
            {
                using (FileStream fs = new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.Read))
                    _mzMLConnection = (Generated.mzMLType)MzmlMethods.mzmlSerializer.Deserialize(fs);
            }

            SourceFile sourceFile;
            if (_mzMLConnection.fileDescription.sourceFileList != null && _mzMLConnection.fileDescription.sourceFileList.sourceFile != null && _mzMLConnection.fileDescription.sourceFileList.sourceFile[0] != null && _mzMLConnection.fileDescription.sourceFileList.sourceFile[0].cvParam != null)
            {
                var simpler = _mzMLConnection.fileDescription.sourceFileList.sourceFile[0];
                string nativeIdFormat = null;
                string fileFormat = null;
                string checkSum = null;
                string checkSumType = null;
                foreach (var cv in simpler.cvParam)
                {
                    if (cv.accession.Equals(@"MS:1000563"))
                    {
                        fileFormat = "Thermo RAW format";
                    }
                    if (cv.accession.Equals(@"MS:1000584"))
                    {
                        fileFormat = "mzML format";
                    }

                    if (cv.accession.Equals(@"MS:1000768"))
                    {
                        nativeIdFormat = "Thermo nativeID format";
                    }
                    if (cv.accession.Equals(@"MS:1000776"))
                    {
                        nativeIdFormat = "scan number only nativeID format";
                    }
                    if (cv.accession.Equals(@"MS:1000824"))
                    {
                        nativeIdFormat = "no nativeID format";
                    }

                    if (cv.accession.Equals(@"MS:1000568"))
                    {
                        checkSum = cv.value;
                        checkSumType = "MD5";
                    }
                    if (cv.accession.Equals(@"MS:1000569"))
                    {
                        checkSum = cv.value;
                        checkSumType = "SHA-1";
                    }
                }

                sourceFile = new SourceFile(
                    nativeIdFormat,
                    fileFormat,
                    checkSum,
                    checkSumType,
                    new Uri(simpler.location),
                    simpler.id,
                    simpler.name);
            }
            else
            {
                string sendCheckSum;
                using (FileStream stream = File.OpenRead(filePath))
                {
                    using (SHA1Managed sha = new SHA1Managed())
                    {
                        byte[] checksum = sha.ComputeHash(stream);
                        sendCheckSum = BitConverter.ToString(checksum)
                            .Replace("-", string.Empty);
                    }
                }
                sourceFile = new SourceFile(
                    @"no nativeID format",
                    @"mzML format",
                    sendCheckSum,
                    @"SHA-1",
                    Path.GetFullPath(filePath),
                    Path.GetFileNameWithoutExtension(filePath));
            }

            var numSpecta = _mzMLConnection.run.spectrumList.spectrum.Length;
            MsDataScan[] scans = new MsDataScan[numSpecta];

            Parallel.ForEach(Partitioner.Create(0, numSpecta), fff =>
            {
                for (int i = fff.Item1; i < fff.Item2; i++)
                {
                    scans[i] = GetMsDataOneBasedScanFromConnection(_mzMLConnection, i + 1, filterParams);
                }
            });

            return new Mzml(scans, sourceFile);
        }

        #endregion Public Methods

    }
}
