using MzLibUtil;

namespace Readers
{
    public class FileReader
    {
        public static TResultFile ReadFile<TResultFile>(string filePath) where TResultFile : IResultFile, new()
        {
            var resultFile = new TResultFile() { FilePath = filePath };
            resultFile.LoadResults();
            return resultFile;
        }

        /// <summary>
        /// Attempts to read in a results file as an IResultFile
        /// </summary>
        /// <param name="filePath">Path to the results file</param>
        /// <returns>An IResultFile created from filePath</returns>
        /// <exception cref="FileNotFoundException">Thrown if filePath doesn't point to a real file</exception>
        /// <exception cref="MzLibException">Thrown if file type is not recognized or can't be converted to IResultFile</exception>
        public static IResultFile ReadResultFile(string filePath)
        {
            if (!File.Exists(filePath) && !Directory.Exists(filePath)) // File and Directory allows Bruker's .d to also work here. 
                throw new FileNotFoundException();
            var resultFileType = filePath.GetResultFileType(); // These calls can throw MzLibExceptions

            // Activator requires an empty constructor, they are guaranteed for any derived class of Readers.ResultFile
            object? resultFile = Activator.CreateInstance(resultFileType);
            if (resultFile is IResultFile castResultFile)
            {
                castResultFile.FilePath = filePath;
                return castResultFile;
            }
            throw new MzLibException($"{resultFileType} files cannot be converted to IResultFile. File path: {filePath}");
        }

        /// <summary>
        /// Attempts to read in a results file as an IQuantifiableResultFile
        /// </summary>
        /// <param name="filePath">Path to the results file</param>
        /// <returns>An IQuantifiableResultFile created from filePath</returns>
        /// <exception cref="FileNotFoundException">Thrown if filePath doesn't point to a real file</exception>
        /// <exception cref="MzLibException">Thrown if file type is not recognized or can't be converted to IQuantifiableResultFile</exception>
        public static IQuantifiableResultFile ReadQuantifiableResultFile(string filePath)
        {
            if (!File.Exists(filePath))
                throw new FileNotFoundException();
            var resultFileType = filePath.GetResultFileType();
            object resultFile = Activator.CreateInstance(resultFileType);
            if (resultFile is IQuantifiableResultFile castResultFile)
            {
                castResultFile.FilePath = filePath;
                return castResultFile;
            }
            throw new MzLibException($"{resultFileType} files cannot be converted to IQuantifiableResultFile");
        }
    }
}
