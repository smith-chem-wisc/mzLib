using MzLibUtil;

namespace Readers
{
    public enum SupportedFileType
    {
        Ms1Feature,
        Ms2Feature,
        Mzrt_TopFd,
        Ms1Tsv_FlashDeconv,
        Tsv_FlashDeconv,
        ThermoRaw,
        MzML,
        Mgf,
        BrukerD,
        psmtsv,
        //osmtsv
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
        CruxResult
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
                SupportedFileType.Mzrt_TopFd => ".mzrt.csv",
                SupportedFileType.Ms1Tsv_FlashDeconv => "_ms1.tsv",
                SupportedFileType.Tsv_FlashDeconv => ".tsv",
                SupportedFileType.ThermoRaw => ".raw",
                SupportedFileType.MzML => ".mzML",
                SupportedFileType.Mgf => ".mgf",
                SupportedFileType.BrukerD => ".d",
                SupportedFileType.psmtsv => ".psmtsv",
                //SupportedFileType.osmtsv => ".osmtsv",
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
                case ".d": return SupportedFileType.BrukerD;
                case ".psmtsv": return SupportedFileType.psmtsv;
                //case ".osmtsv": return SupportedFileType.osmtsv;
                case ".feature":
                    if (filePath.EndsWith(SupportedFileType.Ms1Feature.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.Ms1Feature;
                    if (filePath.EndsWith(SupportedFileType.Ms2Feature.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.Ms2Feature;
                    throw new MzLibException("Feature file type not supported");

                case ".csv":
                    if (filePath.EndsWith(SupportedFileType.Mzrt_TopFd.GetFileExtension(), StringComparison.InvariantCultureIgnoreCase))
                        return SupportedFileType.Mzrt_TopFd;
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

                default:
                    throw new MzLibException("File type not supported");
            }
        }
    }
}
