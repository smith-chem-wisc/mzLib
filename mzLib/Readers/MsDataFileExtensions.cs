using Chemistry;
using MzLibUtil;
using MassSpectrometry; 
namespace Readers
{
    public static class MsDataFileExtensions
    {
        /// <summary>
        /// Export any MsDataFile as an MzML file to a specific file location
        /// CAUTION: some software will check the NativeID for scan numbers
        ///     be sure to set the new NativeID in each MsDataScan if the order has been changed
        /// </summary>
        /// <param name="file"></param>
        /// <param name="destinationPath"></param>
        /// <param name="writeIndexed"></param>
        public static void ExportAsMzML(this MsDataFile file, string destinationPath, bool writeIndexed)
        {
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(file, destinationPath, writeIndexed);
        }
    }
}