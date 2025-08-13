using MzLibUtil;

namespace Readers
{
    public enum SupportedFileType
    {
        Ms1Feature,
        Ms2Feature,
        TopFDMzrt,
        Ms1Tsv_FlashDeconv,
        Tsv_FlashDeconv,
        Tsv_Dinosaur,
        ThermoRaw,
        MzML,
        Mgf,
        Ms1Align,
        Ms2Align,
        psmtsv,
        osmtsv,
        ToppicPrsm,
        ToppicPrsmSingle,
        ToppicProteoform,
        ToppicProteoformSingle,
        MsFraggerPsm,
        MsFraggerPeptide,
        MsFraggerProtein,
        FlashLFQQuantifiedPeak,
        MsPathFinderTTargets,
        MsPathFinderTDecoys,
        MsPathFinderTAllResults,
        CruxResult,
        ExperimentAnnotation,
        BrukerD,
        BrukerTimsTof
    }

    public static class SupportedFileTypeExtensions
    {
        /// <summary>
        /// Returns the extension for the file type
        /// </summary>
        /// <param name="type"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentOutOfRangeException"></exception>
        public static string GetFileExtension(this SupportedFileType type)
        {
            return type switch
            {
                SupportedFileType.Ms1Feature => "_ms1.feature",
                SupportedFileType.Ms2Feature => "_ms2.feature",
                SupportedFileType.TopFDMzrt => ".mzrt.csv",
                SupportedFileType.Ms1Tsv_FlashDeconv => "_ms1.tsv",
                SupportedFileType.Tsv_Dinosaur => ".feature.tsv",
                SupportedFileType.Tsv_FlashDeconv => ".tsv",
                SupportedFileType.ThermoRaw => ".raw",
                SupportedFileType.MzML => ".mzML",
                SupportedFileType.Mgf => ".mgf",
                SupportedFileType.BrukerD => ".d",
                SupportedFileType.BrukerTimsTof => ".d",
                SupportedFileType.Ms1Align => "_ms1.msalign",
                SupportedFileType.Ms2Align => "_ms2.msalign",
                SupportedFileType.psmtsv => ".psmtsv",
                SupportedFileType.osmtsv => ".osmtsv",
                SupportedFileType.ToppicPrsm => "_prsm.tsv",
                SupportedFileType.ToppicPrsmSingle => "_prsm_single.tsv",
                SupportedFileType.ToppicProteoform => "_proteoform.tsv",
                SupportedFileType.ToppicProteoformSingle => "_proteoform_single.tsv",
                SupportedFileType.MsFraggerPsm => "psm.tsv",
                SupportedFileType.MsFraggerPeptide => "peptide.tsv",
                SupportedFileType.MsFraggerProtein => "protein.tsv",
                SupportedFileType.FlashLFQQuantifiedPeak => "Peaks.tsv",
                SupportedFileType.MsPathFinderTTargets => "_IcTarget.tsv",
                SupportedFileType.MsPathFinderTDecoys => "_IcDecoy.tsv",
                SupportedFileType.MsPathFinderTAllResults => "_IcTDA.tsv",
                SupportedFileType.CruxResult => ".txt",
                SupportedFileType.ExperimentAnnotation => "experiment_annotation.tsv",
                _ => throw new MzLibException("File type not supported")
            };
        }
        public static SupportedFileType ParseFileType(this string filePath)
        {
            switch (Path.GetExtension(filePath).ToLower())
            {
                case ".raw": return SupportedFileType.ThermoRaw;
                case ".mzml": return SupportedFileType.MzML;
                case ".mgf": return SupportedFileType.Mgf;
                case ".d":
                    if(!Directory.Exists(filePath)) throw new FileNotFoundException();
                    var fileList = Directory.GetFiles(filePath).Select(p => Path.GetFileName(p));
                    if (fileList.Any(file => file == "analysis.baf"))
                        return SupportedFileType.BrukerD;
                    if (fileList.Any(file => file == "analysis.tdf" || file == "analysis.tsf"))
                        return SupportedFileType.BrukerTimsTof;
                    throw new MzLibException("Bruker file type not recognized");

                case ".psmtsv":
                case ".tsv" when filePath.Contains("Intralinks"):
                    return SupportedFileType.psmtsv;

                case ".osmtsv": 
                    return SupportedFileType.osmtsv;

                case ".feature":
                    if (filePath.EndsWith(SupportedFileType.Ms1Feature.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.Ms1Feature;
                    if (filePath.EndsWith(SupportedFileType.Ms2Feature.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.Ms2Feature;
                    throw new MzLibException("Feature file type not supported");

                case ".csv":
                    if (filePath.EndsWith(SupportedFileType.TopFDMzrt.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.TopFDMzrt;
                    throw new MzLibException("Csv file type not supported");

                case ".tsv":
                {
                    // these tsv cases have a specialized ending before the .tsv
                    if (filePath.EndsWith(SupportedFileType.Ms1Tsv_FlashDeconv.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.Ms1Tsv_FlashDeconv;
                    if (filePath.EndsWith(SupportedFileType.ToppicPrsm.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.ToppicPrsm;
                    if (filePath.EndsWith(SupportedFileType.ToppicProteoform.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.ToppicProteoform;
                    if (filePath.EndsWith(SupportedFileType.ToppicPrsmSingle.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.ToppicPrsmSingle;
                    if (filePath.EndsWith(SupportedFileType.ToppicProteoformSingle.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.ToppicProteoformSingle;
                    if (filePath.EndsWith(SupportedFileType.MsFraggerPsm.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.MsFraggerPsm;
                    if (filePath.EndsWith(SupportedFileType.MsFraggerPeptide.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.MsFraggerPeptide;
                    if (filePath.EndsWith(SupportedFileType.MsFraggerProtein.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.MsFraggerProtein;
                    if (filePath.EndsWith(SupportedFileType.FlashLFQQuantifiedPeak.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.FlashLFQQuantifiedPeak;
                    if (filePath.EndsWith(SupportedFileType.MsPathFinderTTargets.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.MsPathFinderTTargets;
                    if (filePath.EndsWith(SupportedFileType.MsPathFinderTDecoys.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.MsPathFinderTDecoys;
                    if (filePath.EndsWith(SupportedFileType.MsPathFinderTAllResults.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.MsPathFinderTAllResults;
                    if(filePath.EndsWith(SupportedFileType.ExperimentAnnotation.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.ExperimentAnnotation;
                    if(filePath.EndsWith(SupportedFileType.Tsv_Dinosaur.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.Tsv_Dinosaur;

                        // these tsv cases are just .tsv and need an extra step to determine the type
                        // currently need to distinguish between FlashDeconvTsv and MsFraggerPsm
                        using var sw = new StreamReader(filePath);
                    var firstLine = sw.ReadLine() ?? "";
                    if (firstLine == "") throw new MzLibException("Tsv file is empty");

                    if (firstLine.Contains("FeatureIndex"))
                        return SupportedFileType.Tsv_FlashDeconv;
                    throw new MzLibException("Tsv file type not supported");
                }

                case ".txt":
                    if (filePath.EndsWith(SupportedFileType.CruxResult.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.CruxResult;
                    throw new MzLibException("Txt file type not supported");

                case ".msalign":
                    if (filePath.EndsWith(SupportedFileType.Ms1Align.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.Ms1Align;
                    if (filePath.EndsWith(SupportedFileType.Ms2Align.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.Ms2Align;
                    throw new MzLibException("MsAlign file type not supported, must end with _msX.msalign where X is 1 or 2");

                default:
                    throw new MzLibException("File type not supported");
            }
        }

        /// <summary>
        /// Returns the typeOf the related class by parsing the SupportedFileType enum
        /// If the SupportedFileType is not an IResultFile, an MzLibException is thrown
        /// </summary>
        /// <param name="type"></param>
        /// <returns>Type of File Reader class</returns>
        /// <exception cref="MzLibException">Throws exception if SupportedFileType is unrecognized or is not an IResultFile</exception>
        public static Type GetResultFileType(this SupportedFileType type)
        {
            return type switch
            {
                SupportedFileType.Ms1Feature => typeof(Ms1FeatureFile),
                SupportedFileType.Ms2Feature => typeof(Ms2FeatureFile),
                SupportedFileType.TopFDMzrt => typeof(TopFDMzrtFile),
                SupportedFileType.Ms1Tsv_FlashDeconv => typeof(FlashDeconvMs1TsvFile),
                SupportedFileType.Tsv_FlashDeconv => typeof(FlashDeconvTsvFile),
                SupportedFileType.Tsv_Dinosaur => typeof(DinosaurTsvFile),
                SupportedFileType.psmtsv => typeof(PsmFromTsvFile),
                SupportedFileType.osmtsv => typeof(OsmFromTsvFile),
                SupportedFileType.ToppicPrsm => typeof(ToppicSearchResultFile),
                SupportedFileType.ToppicPrsmSingle => typeof(ToppicSearchResultFile),
                SupportedFileType.ToppicProteoform => typeof(ToppicSearchResultFile),
                SupportedFileType.ToppicProteoformSingle => typeof(ToppicSearchResultFile),
                SupportedFileType.MsFraggerPsm => typeof(MsFraggerPsmFile),
                SupportedFileType.MsFraggerPeptide => typeof(MsFraggerPeptideFile),
                SupportedFileType.MsFraggerProtein => typeof(MsFraggerProteinFile),
                SupportedFileType.FlashLFQQuantifiedPeak => typeof(QuantifiedPeakFile),
                SupportedFileType.MsPathFinderTTargets => typeof(MsPathFinderTResultFile),
                SupportedFileType.MsPathFinderTDecoys => typeof(MsPathFinderTResultFile),
                SupportedFileType.MsPathFinderTAllResults => typeof(MsPathFinderTResultFile),
                SupportedFileType.CruxResult => typeof(CruxResultFile),
                SupportedFileType.ExperimentAnnotation => typeof(ExperimentAnnotationFile),
                SupportedFileType.ThermoRaw => typeof(MsDataFileToResultFileAdapter),
                SupportedFileType.MzML => typeof(MsDataFileToResultFileAdapter),
                SupportedFileType.Mgf => typeof(MsDataFileToResultFileAdapter),
                SupportedFileType.BrukerD => typeof(MsDataFileToResultFileAdapter),
                SupportedFileType.BrukerTimsTof => typeof(MsDataFileToResultFileAdapter),
                SupportedFileType.Ms1Align => typeof(MsDataFileToResultFileAdapter),
                SupportedFileType.Ms2Align => typeof(MsDataFileToResultFileAdapter),
                _ => throw new MzLibException("File type not supported")
            };
        }

        public static Type GetResultFileType(this string filePath) => filePath.ParseFileType().GetResultFileType();
    }
}
