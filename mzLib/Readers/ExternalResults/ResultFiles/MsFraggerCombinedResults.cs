using CsvHelper;
using Readers.ExternalResults.BaseClasses;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using MathNet.Numerics;

namespace Readers
{
    public class MsFraggerCombinedResults : ResultFile<MsFraggerPsm>, IResultFile, IQuantifiableResultFile
    {
        #region Properties/Fields

        public string FullFolderPath => FilePath;
        private List<string> allPsmFilePaths;
        public List<MsFraggerPsmFile> AllPsmFiles { get; private set; }
        public ExperimentAnnotationFile ExperimentAnnotations { get; private set; }

        #endregion

        #region IResultFile Implementatation

        public override SupportedFileType FileType => SupportedFileType.MsFraggerPsm;
        public override Software Software { get; set; }
        public MsFraggerCombinedResults(string filePath) : base(filePath, Software.MsFragger) { }

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
            throw new Exception("Method not yet implemented.");
        }

        #endregion

        public void LoadExperimentAnnotationResults()
        {
            string combinedFilePath = Path.Combine(FullFolderPath, "experiment_annotation.tsv");
            if (!File.Exists(combinedFilePath)) { throw new FileNotFoundException("The experiment_annotation.tsv file was not found"); }

            ExperimentAnnotations = new ExperimentAnnotationFile(combinedFilePath);
        }

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

        public Dictionary<string, string> FileNameToFilePath(List<string> filePaths)
        {
            filePaths = filePaths.Distinct().ToList();
            List<string> fileNames = Results.Select(psm => psm.FileName).Distinct().ToList();
            Dictionary<string, string> allFiles = new Dictionary<string, string>();

            foreach (var name in fileNames)
            {
                string fileName = Path.GetFileNameWithoutExtension(name);
                if (fileName.Contains("."))
                {
                    fileName = Path.GetFileNameWithoutExtension(fileName);
                }

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

        public Dictionary<string, string> FileNameToFilePath()
        {
            List<string> filePaths = ExperimentAnnotations.Select(psm => psm.File).Distinct().ToList();
            List<string> fileNames = Results.Select(psm => psm.FileName).Distinct().ToList();
            Dictionary<string, string> allFiles = new Dictionary<string, string>();

            foreach (var name in fileNames)
            {
                string fileName = Path.GetFileNameWithoutExtension(name);
                if (fileName.Contains("."))
                {
                    fileName = Path.GetFileNameWithoutExtension(fileName);
                }

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
