namespace Readers
{
    public class MsFraggerCombinedResults : ResultFile<MsFraggerPsm>, IQuantifiableResultFile
    {
        #region Properties/Fields

        public string FullFolderPath => FilePath; // The full file path to the folder of MSFragger results
        private List<string> allPsmFilePaths; // List of the full file paths to the psm files of every sample

        // A list of all the MSFraggerPsmFile objects that correspond to each sample within an experiment
        public List<MsFraggerPsmFile> AllPsmFiles { get; private set; }

        // Contains descriptive information on every ms data file in the experiment (sample name, full path to the ms data file, etc.)
        public ExperimentAnnotationFile ExperimentAnnotations { get; private set; }

        #endregion

        #region IResultFile Implementatation

        public override SupportedFileType FileType => SupportedFileType.MsFraggerPsm;
        public override Software Software { get; set; }
        public MsFraggerCombinedResults(string filePath) : base(filePath, Software.MsFragger) { }

        /// <summary>
        /// Loads the results from each psm.tsv file in the results folder, builds one list of MsFraggerPsms,
        /// and Calls LoadExperimentAnnotation, FindAllFilePaths, LoadPsmResults,
        /// then selects every results from each MsFraggerPsmFile in AllPsmFiles and writes them to one concatenated list. 
        /// </summary>
        public override void LoadResults()
        {
            LoadExperimentAnnotationResults();
            FindAllFilePaths();
            LoadPsmResults();

            List<MsFraggerPsm> concatList = new List<MsFraggerPsm>();
            foreach (var file in AllPsmFiles)
            {
                concatList.AddRange(file);
            }

            Results = concatList;
        }

        public override void WriteResults(string outputPath)
        {
            throw new NotImplementedException("Method not yet implemented.");
        }

        #endregion

        /// <summary>
        /// Checks for existence of experiment annotation file and loads its it as an ExperimentAnnotationResultFile, 
        /// then sets the ExperimentAnnotations property
        /// </summary>
        /// <exception cref="FileNotFoundException"></exception>
        public void LoadExperimentAnnotationResults()
        {
            string combinedFilePath = Path.Combine(FullFolderPath, "experiment_annotation.tsv");
            if (!File.Exists(combinedFilePath)) { throw new FileNotFoundException("The experiment_annotation.tsv file was not found"); }

            ExperimentAnnotations = new ExperimentAnnotationFile(combinedFilePath);
        }

        /// <summary>
        /// For each path in AllPsmFilePaths, creates and loads an MsFraggerPsmFile. 
        /// Then constructs the AllPsmFiles list
        /// </summary>
        public void LoadPsmResults()
        {
            AllPsmFiles = new List<MsFraggerPsmFile>();

            foreach(var path in allPsmFilePaths)
            {
                MsFraggerPsmFile file = new MsFraggerPsmFile(path);
                AllPsmFiles.Add(file);
            }
        }

        public IEnumerable<IQuantifiableRecord> GetQuantifiableResults() => Results;

        /// <summary>
        /// Links the file name associated with the an IQuantifiableRecord 
        /// to the raw file path of MassSpec data in the fullFilePath list
        /// </summary>
        /// <param name="filePaths"> list of file paths associated with each distinct record </param>
        /// <returns> Dictionary of file names and their associted full paths </returns>
        public Dictionary<string, string> FileNameToFilePath(List<string> filePaths)
        {
            Dictionary<string, string> allFiles = new Dictionary<string, string>();

            allFiles = AllPsmFiles.Select(file => file.FileNameToFilePath(filePaths))
                                  .SelectMany(dictionary => dictionary)
                                  .GroupBy(x => x.Key)
                                  .Select(keyValuePair => keyValuePair.First())
                                  .ToDictionary(fileName => fileName.Key, filePath => filePath.Value);

            return allFiles;
        }

        /// <summary>
        /// Links the file name associated with IQuantifiableRecord to the raw file path pf MassSpec file 
        /// using the full file paths from the experiment annotation file.
        /// </summary>
        /// <returns> Dictionary of file names and their associted full paths </returns>
        public Dictionary<string, string> FileNameToFilePath()
        {
            List<string> filePaths = ExperimentAnnotations.Select(psm => psm.File).Distinct().ToList();
            List<string> fileNames = Results.Select(psm => psm.FileName).Distinct().ToList();
            Dictionary<string, string> allFiles = new Dictionary<string, string>();

            foreach (var name in fileNames)
            {
                string fileName = Path.GetFileName(name);

                // MSFragger results append the raw file with "interact-" and replace .raw with .pep.xml
                // In order to correctly match the file names, these changes must be removed
                fileName = fileName.Replace("interact-", "").Replace(".pep.xml", "");

                foreach (var path in filePaths)
                {
                    if (path.Contains(fileName) && !allFiles.ContainsKey(name))
                    {
                        allFiles.Add(name, path);
                        break;
                    }
                }
            }

            return allFiles;
        }

        /// <summary>
        /// Uses the ExperimentAnnotations to locate each psm.tsv file in the results folder. 
        /// Adds the path to each psm.tsv file in the results folder to AllPsmFilePaths
        /// </summary>
        /// <exception cref="FileNotFoundException"></exception>
        private void FindAllFilePaths()
        {
            allPsmFilePaths = new List<string>();

            List<string> sampleNames = ExperimentAnnotations.Select(psm => psm.SampleName).Distinct().ToList();
            string[] directoryEntries = Directory.GetDirectories(FullFolderPath);

            foreach (var directoryEntry in directoryEntries)
            {
                string directoryName = Path.GetFileName(directoryEntry.TrimEnd(Path.DirectorySeparatorChar));

                foreach (var sample in sampleNames)
                {
                    if (directoryName.Equals(sample))
                    {
                        string psmFile = Path.Combine(directoryEntry, "psm.tsv");
                        if (!File.Exists(psmFile)) { throw new FileNotFoundException("This psm.tsv file was not found"); }

                        allPsmFilePaths.Add(psmFile);
                    }
                }
            }
        }
    }
}