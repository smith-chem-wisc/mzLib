﻿using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Xml.Serialization;

namespace IO.MzML
{
    public static class MzmlMethods
    {

        #region Internal Fields

        internal static readonly XmlSerializer indexedSerializer = new XmlSerializer(typeof(Generated.indexedmzML));
        internal static readonly XmlSerializer mzmlSerializer = new XmlSerializer(typeof(Generated.mzMLType));

        #endregion Internal Fields

        #region Private Fields

        private static readonly Dictionary<DissociationType, string> DissociationTypeAccessions = new Dictionary<DissociationType, string>{
            {DissociationType.CID, "MS:1000133"},
            {DissociationType.ISCID, "MS:1001880"},
            {DissociationType.HCD, "MS:1000422" },
            {DissociationType.ETD, "MS:1000598"},
            {DissociationType.MPD, "MS:1000435"},
            {DissociationType.PQD, "MS:1000599"},
            {DissociationType.Unknown, "MS:1000044"} };

        private static readonly Dictionary<DissociationType, string> DissociationTypeNames = new Dictionary<DissociationType, string>{
            {DissociationType.CID, "collision-induced dissociation"},
            {DissociationType.ISCID, "in-source collision-induced dissociation"},
            {DissociationType.HCD, "beam-type collision-induced dissociation"},
            {DissociationType.ETD, "electron transfer dissociation"},
            {DissociationType.MPD, "photodissociation"},
            {DissociationType.PQD, "pulsed q dissociation"},
            {DissociationType.Unknown, "dissociation method"}};

        private static readonly Dictionary<bool, string> CentroidAccessions = new Dictionary<bool, string>{
            {true, "MS:1000127"},
            {false, "MS:1000128"}};

        private static readonly Dictionary<bool, string> CentroidNames = new Dictionary<bool, string>{
            {true, "centroid spectrum"},
            {false, "profile spectrum"}};

        private static readonly Dictionary<Polarity, string> PolarityAccessions = new Dictionary<Polarity, string>{
            {Polarity.Negative, "MS:1000129"},
            {Polarity.Positive, "MS:1000130"}};

        private static readonly Dictionary<Polarity, string> PolarityNames = new Dictionary<Polarity, string>{
            {Polarity.Negative, "negative scan"},
            {Polarity.Positive, "positive scan"}};
        #endregion Private Fields

        #region Public Methods

        public static void CreateAndWriteMyMzmlWithCalibratedSpectra(IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile, string outputFile, bool writeIndexed)
        {
            var mzML = new Generated.mzMLType()
            {
                version = "1",
                cvList = new Generated.CVListType()
            };
            mzML.cvList.count = "1";
            mzML.cvList.cv = new Generated.CVType[1];
            mzML.cvList.cv[0] = new Generated.CVType()
            {
                URI = @"https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo",
                fullName = "Proteomics Standards Initiative Mass Spectrometry Ontology",
                id = "MS"
            };
            mzML.fileDescription = new Generated.FileDescriptionType()
            {
                fileContent = new Generated.ParamGroupType()
            };
            mzML.fileDescription.fileContent.cvParam = new Generated.CVParamType[2];
            mzML.fileDescription.fileContent.cvParam[0] = new Generated.CVParamType()
            {
                accession = "MS:1000579" // MS1 Data
            };
            mzML.fileDescription.fileContent.cvParam[1] = new Generated.CVParamType()
            {
                accession = "MS:1000580" // MSn Data
            };
            mzML.softwareList = new Generated.SoftwareListType()
            {
                count = "1",
                software = new Generated.SoftwareType[1]
            };

            // TODO: add the raw file fields
            mzML.softwareList.software[0] = new Generated.SoftwareType()
            {
                id = "mzLib",
                version = "1",
                cvParam = new Generated.CVParamType[1]
            };
            mzML.softwareList.software[0].cvParam[0] = new Generated.CVParamType()
            {
                accession = "MS:1000799",
                value = "mzLib"
            };

            // Leaving empty. Can't figure out the configurations.
            // ToDo: read instrumentConfigurationList from mzML file
            mzML.instrumentConfigurationList = new Generated.InstrumentConfigurationListType();

            mzML.dataProcessingList = new Generated.DataProcessingListType()
            {
                count = "1",
                dataProcessing = new Generated.DataProcessingType[1]
            };
            // Only writing mine! Might have had some other data processing (but not if it is a raw file)
            // ToDo: read dataProcessingList from mzML file
            mzML.dataProcessingList.dataProcessing[0] = new Generated.DataProcessingType()
            {
                id = "mzLibProcessing"
            };
            mzML.run = new Generated.RunType()
            {
                chromatogramList = new Generated.ChromatogramListType()
                {
                    count = "1",
                    chromatogram = new Generated.ChromatogramType[1]
                }
            };
            // ToDo: Finish the chromatogram writing!
            mzML.run.chromatogramList.chromatogram[0] = new Generated.ChromatogramType();

            mzML.run.spectrumList = new Generated.SpectrumListType()
            {
                count = (myMsDataFile.NumSpectra).ToString(CultureInfo.InvariantCulture),
                defaultDataProcessingRef = "mzLibProcessing",
                spectrum = new Generated.SpectrumType[myMsDataFile.NumSpectra]
            };

            // Loop over all spectra
            for (int i = 1; i <= myMsDataFile.NumSpectra; i++)
            {
                mzML.run.spectrumList.spectrum[i - 1] = new Generated.SpectrumType()
                {
                    defaultArrayLength = myMsDataFile.GetOneBasedScan(i).MassSpectrum.Size,

                    index = i.ToString(CultureInfo.InvariantCulture),
                    id = myMsDataFile.GetOneBasedScan(i).OneBasedScanNumber.ToString(),

                    cvParam = new Generated.CVParamType[8]
                };
                mzML.run.spectrumList.spectrum[i - 1].cvParam[0] = new Generated.CVParamType();

                if (myMsDataFile.GetOneBasedScan(i).MsnOrder == 1)
                {
                    mzML.run.spectrumList.spectrum[i - 1].cvParam[0].accession = "MS:1000579";
                }
                else if (myMsDataFile.GetOneBasedScan(i) is IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>)
                {
                    var scanWithPrecursor = myMsDataFile.GetOneBasedScan(i) as IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>;
                    mzML.run.spectrumList.spectrum[i - 1].cvParam[0].accession = "MS:1000580";

                    // So needs a precursor!
                    mzML.run.spectrumList.spectrum[i - 1].precursorList = new Generated.PrecursorListType()
                    {
                        count = 1.ToString(),
                        precursor = new Generated.PrecursorType[1]
                    };
                    mzML.run.spectrumList.spectrum[i - 1].precursorList.precursor[0] = new Generated.PrecursorType();
                    string precursorID = scanWithPrecursor.OneBasedPrecursorScanNumber.ToString();
                    mzML.run.spectrumList.spectrum[i - 1].precursorList.precursor[0].spectrumRef = precursorID;
                    mzML.run.spectrumList.spectrum[i - 1].precursorList.precursor[0].selectedIonList = new Generated.SelectedIonListType()
                    {
                        count = 1.ToString(),
                        selectedIon = new Generated.ParamGroupType[1]
                    };
                    mzML.run.spectrumList.spectrum[i - 1].precursorList.precursor[0].selectedIonList.selectedIon[0] = new Generated.ParamGroupType()
                    {
                        cvParam = new Generated.CVParamType[3]
                    };
                    // Selected ion MZ
                    mzML.run.spectrumList.spectrum[i - 1].precursorList.precursor[0].selectedIonList.selectedIon[0].cvParam[0] = new Generated.CVParamType()
                    {
                        name = "selected ion m/z",
                        value = scanWithPrecursor.SelectedIonMZ.ToString(CultureInfo.InvariantCulture),
                        accession = "MS:1000744"
                    };

                    // Charge State
                    if (scanWithPrecursor.SelectedIonChargeStateGuess.HasValue)
                    {
                        mzML.run.spectrumList.spectrum[i - 1].precursorList.precursor[0].selectedIonList.selectedIon[0].cvParam[1] = new Generated.CVParamType()
                        {
                            name = "charge state",
                            value = scanWithPrecursor.SelectedIonChargeStateGuess.Value.ToString(CultureInfo.InvariantCulture),
                            accession = "MS:1000041"
                        };
                    }

                    // Selected ion intensity
                    if (scanWithPrecursor.SelectedIonIntensity.HasValue)
                    {
                        mzML.run.spectrumList.spectrum[i - 1].precursorList.precursor[0].selectedIonList.selectedIon[0].cvParam[2] = new Generated.CVParamType()
                        {
                            name = "peak intensity",
                            value = scanWithPrecursor.SelectedIonIntensity.Value.ToString(CultureInfo.InvariantCulture),
                            accession = "MS:1000042"
                        };
                    }

                    MzRange isolationRange = scanWithPrecursor.IsolationRange;
                    mzML.run.spectrumList.spectrum[i - 1].precursorList.precursor[0].isolationWindow = new Generated.ParamGroupType()
                    {
                        cvParam = new Generated.CVParamType[3]
                    };
                    mzML.run.spectrumList.spectrum[i - 1].precursorList.precursor[0].isolationWindow.cvParam[0] = new Generated.CVParamType()
                    {
                        accession = "MS:1000827",
                        name = "isolation window target m/z",
                        value = isolationRange.Mean.ToString(CultureInfo.InvariantCulture)
                    };
                    mzML.run.spectrumList.spectrum[i - 1].precursorList.precursor[0].isolationWindow.cvParam[1] = new Generated.CVParamType()
                    {
                        accession = "MS:1000828",
                        name = "isolation window lower offset",
                        value = (isolationRange.Width / 2).ToString(CultureInfo.InvariantCulture)
                    };
                    mzML.run.spectrumList.spectrum[i - 1].precursorList.precursor[0].isolationWindow.cvParam[2] = new Generated.CVParamType()
                    {
                        accession = "MS:1000829",
                        name = "isolation window upper offset",
                        value = (isolationRange.Width / 2).ToString(CultureInfo.InvariantCulture)
                    };
                    mzML.run.spectrumList.spectrum[i - 1].precursorList.precursor[0].activation = new Generated.ParamGroupType()
                    {
                        cvParam = new Generated.CVParamType[1]
                    };
                    mzML.run.spectrumList.spectrum[i - 1].precursorList.precursor[0].activation.cvParam[0] = new Generated.CVParamType();

                    DissociationType dissociationType = scanWithPrecursor.DissociationType;

                    mzML.run.spectrumList.spectrum[i - 1].precursorList.precursor[0].activation.cvParam[0].accession = DissociationTypeAccessions[dissociationType];
                    mzML.run.spectrumList.spectrum[i - 1].precursorList.precursor[0].activation.cvParam[0].name = DissociationTypeNames[dissociationType];
                }

                mzML.run.spectrumList.spectrum[i - 1].cvParam[1] = new Generated.CVParamType()
                {
                    name = "ms level",
                    accession = "MS:1000511",
                    value = myMsDataFile.GetOneBasedScan(i).MsnOrder.ToString(CultureInfo.InvariantCulture)
                };
                mzML.run.spectrumList.spectrum[i - 1].cvParam[2] = new Generated.CVParamType()
                {
                    name = CentroidNames[myMsDataFile.GetOneBasedScan(i).IsCentroid],
                    accession = CentroidAccessions[myMsDataFile.GetOneBasedScan(i).IsCentroid]
                };
                if (PolarityNames.TryGetValue(myMsDataFile.GetOneBasedScan(i).Polarity, out string polarityName) && PolarityAccessions.TryGetValue(myMsDataFile.GetOneBasedScan(i).Polarity, out string polarityAccession))
                {
                    mzML.run.spectrumList.spectrum[i - 1].cvParam[3] = new Generated.CVParamType()
                    {
                        name = polarityName,
                        accession = polarityAccession
                    };
                }
                // Spectrum title
                mzML.run.spectrumList.spectrum[i - 1].cvParam[4] = new Generated.CVParamType()
                {
                    name = "spectrum title",
                    accession = "MS:1000796",
                    value = myMsDataFile.GetOneBasedScan(i).OneBasedScanNumber.ToString()
                };
                if ((myMsDataFile.GetOneBasedScan(i).MassSpectrum.Size) > 0)
                {
                    // Lowest observed mz
                    mzML.run.spectrumList.spectrum[i - 1].cvParam[5] = new Generated.CVParamType()
                    {
                        name = "lowest observed m/z",
                        accession = "MS:1000528",
                        value = myMsDataFile.GetOneBasedScan(i).MassSpectrum.FirstX.ToString(CultureInfo.InvariantCulture)
                    };

                    // Highest observed mz
                    mzML.run.spectrumList.spectrum[i - 1].cvParam[6] = new Generated.CVParamType()
                    {
                        name = "highest observed m/z",
                        accession = "MS:1000527",
                        value = myMsDataFile.GetOneBasedScan(i).MassSpectrum.LastX.ToString(CultureInfo.InvariantCulture)
                    };
                }

                // Total ion current
                mzML.run.spectrumList.spectrum[i - 1].cvParam[7] = new Generated.CVParamType()
                {
                    name = "total ion current",
                    accession = "MS:1000285",
                    value = myMsDataFile.GetOneBasedScan(i).TotalIonCurrent.ToString(CultureInfo.InvariantCulture)
                };

                // Retention time
                mzML.run.spectrumList.spectrum[i - 1].scanList = new Generated.ScanListType()
                {
                    count = "1",
                    scan = new Generated.ScanType[1]
                };
                mzML.run.spectrumList.spectrum[i - 1].scanList.scan[0] = new Generated.ScanType()
                {
                    cvParam = new Generated.CVParamType[3]
                };
                mzML.run.spectrumList.spectrum[i - 1].scanList.scan[0].cvParam[0] = new Generated.CVParamType()
                {
                    name = "scan start time",
                    accession = "MS:1000016",
                    value = myMsDataFile.GetOneBasedScan(i).RetentionTime.ToString(CultureInfo.InvariantCulture),
                    unitCvRef = "UO",
                    unitAccession = "UO:0000031",
                    unitName = "minute"
                };
                mzML.run.spectrumList.spectrum[i - 1].scanList.scan[0].cvParam[1] = new Generated.CVParamType()
                {
                    name = "filter string",
                    accession = "MS:1000512",
                    value = myMsDataFile.GetOneBasedScan(i).ScanFilter
                };
                if (myMsDataFile.GetOneBasedScan(i).InjectionTime.HasValue)
                {
                    mzML.run.spectrumList.spectrum[i - 1].scanList.scan[0].cvParam[2] = new Generated.CVParamType()
                    {
                        name = "ion injection time",
                        accession = "MS:1000927",
                        value = myMsDataFile.GetOneBasedScan(i).InjectionTime.Value.ToString(CultureInfo.InvariantCulture)
                    };
                }
                if (myMsDataFile.GetOneBasedScan(i) is IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>)
                {
                    var scanWithPrecursor = myMsDataFile.GetOneBasedScan(i) as IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>;
                    if (scanWithPrecursor.SelectedIonMonoisotopicGuessMz.HasValue)
                    {
                        mzML.run.spectrumList.spectrum[i - 1].scanList.scan[0].userParam = new Generated.UserParamType[1];
                        mzML.run.spectrumList.spectrum[i - 1].scanList.scan[0].userParam[0] = new Generated.UserParamType()
                        {
                            name = "[mzLib]Monoisotopic M/Z:",
                            value = scanWithPrecursor.SelectedIonMonoisotopicGuessMz.Value.ToString(CultureInfo.InvariantCulture)
                        };
                    }
                }

                mzML.run.spectrumList.spectrum[i - 1].scanList.scan[0].scanWindowList = new Generated.ScanWindowListType()
                {
                    count = 1,
                    scanWindow = new Generated.ParamGroupType[1]
                };
                mzML.run.spectrumList.spectrum[i - 1].scanList.scan[0].scanWindowList.scanWindow[0] = new Generated.ParamGroupType()
                {
                    cvParam = new Generated.CVParamType[2]
                };
                mzML.run.spectrumList.spectrum[i - 1].scanList.scan[0].scanWindowList.scanWindow[0].cvParam[0] = new Generated.CVParamType()
                {
                    name = "scan window lower limit",
                    accession = "MS:1000501",
                    value = myMsDataFile.GetOneBasedScan(i).ScanWindowRange.Minimum.ToString(CultureInfo.InvariantCulture)
                };
                mzML.run.spectrumList.spectrum[i - 1].scanList.scan[0].scanWindowList.scanWindow[0].cvParam[1] = new Generated.CVParamType()
                {
                    name = "scan window upper limit",
                    accession = "MS:1000500",
                    value = myMsDataFile.GetOneBasedScan(i).ScanWindowRange.Maximum.ToString(CultureInfo.InvariantCulture)
                };
                mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList = new Generated.BinaryDataArrayListType()
                {
                    // ONLY WRITING M/Z AND INTENSITY DATA, NOT THE CHARGE! (but can add charge info later)
                    // CHARGE (and other stuff) CAN BE IMPORTANT IN ML APPLICATIONS!!!!!
                    count = 2.ToString(),
                    binaryDataArray = new Generated.BinaryDataArrayType[5]
                };

                // M/Z Data
                mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[0] = new Generated.BinaryDataArrayType()
                {
                    binary = myMsDataFile.GetOneBasedScan(i).MassSpectrum.Get64BitXarray()
                };
                mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[0].encodedLength = (4 * Math.Ceiling(((double)mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[0].binary.Length / 3))).ToString(CultureInfo.InvariantCulture);
                mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[0].cvParam = new Generated.CVParamType[3];
                mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[0].cvParam[0] = new Generated.CVParamType()
                {
                    accession = "MS:1000514",
                    name = "m/z array"
                };
                mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[0].cvParam[1] = new Generated.CVParamType()
                {
                    accession = "MS:1000523",
                    name = "64-bit float"
                };
                mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[0].cvParam[2] = new Generated.CVParamType()
                {
                    accession = "MS:1000576",
                    name = "no compression"
                };

                // Intensity Data
                mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[1] = new Generated.BinaryDataArrayType()
                {
                    binary = myMsDataFile.GetOneBasedScan(i).MassSpectrum.Get64BitYarray()
                };
                mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[1].encodedLength = (4 * Math.Ceiling(((double)mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[1].binary.Length / 3))).ToString(CultureInfo.InvariantCulture);
                mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[1].cvParam = new Generated.CVParamType[3];
                mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[1].cvParam[0] = new Generated.CVParamType()
                {
                    accession = "MS:1000515",
                    name = "intensity array"
                };
                mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[1].cvParam[1] = new Generated.CVParamType()
                {
                    accession = "MS:1000523",
                    name = "64-bit float"
                };
                mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[1].cvParam[2] = new Generated.CVParamType()
                {
                    accession = "MS:1000576",
                    name = "no compression"
                };
                if (myMsDataFile.GetOneBasedScan(i).NoiseData != null)
                {
                    // mass
                    mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[2] = new Generated.BinaryDataArrayType()
                    {
                        binary = myMsDataFile.GetOneBasedScan(i).Get64BitNoiseDataMass()
                    };
                    mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[2].encodedLength = (4 * Math.Ceiling(((double)mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[2].binary.Length / 3))).ToString(CultureInfo.InvariantCulture);
                    mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[2].cvParam = new Generated.CVParamType[3];
                    mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[2].cvParam[0] = new Generated.CVParamType()
                    {
                        accession = "MS:1000786",
                        name = "non-standard data array"
                    };
                    mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[2].cvParam[1] = new Generated.CVParamType()
                    {
                        accession = "MS:1000523",
                        name = "64-bit float"
                    };
                    mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[2].cvParam[2] = new Generated.CVParamType()
                    {
                        accession = "MS:1000576",
                        name = "no compression"
                    };
                    mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[2].userParam = new Generated.UserParamType[1];
                    mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[2].userParam[0] = new Generated.UserParamType()
                    {
                        name = "kelleherCustomType",
                        value = "noise m/z"
                    };

                    // noise
                    mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[3] = new Generated.BinaryDataArrayType()
                    {
                        binary = myMsDataFile.GetOneBasedScan(i).Get64BitNoiseDataNoise()
                    };
                    mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[3].encodedLength = (4 * Math.Ceiling(((double)mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[3].binary.Length / 3))).ToString(CultureInfo.InvariantCulture);
                    mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[3].cvParam = new Generated.CVParamType[3];
                    mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[3].cvParam[0] = new Generated.CVParamType()
                    {
                        accession = "MS:1000786",
                        name = "non-standard data array"
                    };
                    mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[3].cvParam[1] = new Generated.CVParamType()
                    {
                        accession = "MS:1000523",
                        name = "64-bit float"
                    };
                    mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[3].cvParam[2] = new Generated.CVParamType()
                    {
                        accession = "MS:1000576",
                        name = "no compression"
                    };
                    mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[3].userParam = new Generated.UserParamType[1];
                    mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[3].userParam[0] = new Generated.UserParamType()
                    {
                        name = "kelleherCustomType",
                        value = "noise baseline"
                    };

                    // baseline
                    mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[4] = new Generated.BinaryDataArrayType()
                    {
                        binary = myMsDataFile.GetOneBasedScan(i).Get64BitNoiseDataBaseline()
                    };
                    mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[4].encodedLength = (4 * Math.Ceiling(((double)mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[4].binary.Length / 3))).ToString(CultureInfo.InvariantCulture);
                    mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[4].cvParam = new Generated.CVParamType[3];
                    mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[4].cvParam[0] = new Generated.CVParamType()
                    {
                        accession = "MS:1000786",
                        name = "non-standard data array"
                    };
                    mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[4].cvParam[1] = new Generated.CVParamType()
                    {
                        accession = "MS:1000523",
                        name = "64-bit float"
                    };
                    mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[4].cvParam[2] = new Generated.CVParamType()
                    {
                        accession = "MS:1000576",
                        name = "no compression"
                    };
                    mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[4].userParam = new Generated.UserParamType[1];
                    mzML.run.spectrumList.spectrum[i - 1].binaryDataArrayList.binaryDataArray[4].userParam[0] = new Generated.UserParamType()
                    {
                        name = "kelleherCustomType",
                        value = "noise intensity"
                    };
                }
            }

            if (writeIndexed)
                throw new NotImplementedException("Writing indexed mzMLs not yet supported");
            else
            {
                using (TextWriter writer = new StreamWriter(outputFile))
                {
                    mzmlSerializer.Serialize(writer, mzML);
                }
            }
        }

        #endregion Public Methods

    }
}