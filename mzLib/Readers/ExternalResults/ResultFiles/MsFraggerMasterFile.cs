using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers
{
    public class MsFraggerMasterFile
    {
        public MsFraggerExperimentFile ExperimentFile { get; set; }
        public MsFraggerPeptideFile PeptideFile { get; set; }
        public Dictionary<MsFraggerExperiment, MsFraggerPsmFile> PsmFiles { get; set; }

        public static MsFraggerMasterFile Load(string pathToResultsFolder)
        {
            string experimentFilePath = Path.Combine(pathToResultsFolder, "experiment_annotation.tsv");
            string peptideFilePath = Path.Combine(pathToResultsFolder, "combined_modified_peptide.tsv");
            if(!File.Exists(peptideFilePath))
            {
                peptideFilePath = Path.Combine(pathToResultsFolder, "combined_modified_peptide.tsv");
                if(!File.Exists(peptideFilePath))
                {
                    peptideFilePath = null;
                }
            }

            if(!File.Exists(experimentFilePath))
            {
                throw new FileNotFoundException("Could not find the experiment file at " + experimentFilePath);
            }

            var experimentFile = new MsFraggerExperimentFile(experimentFilePath);
            experimentFile.LoadResults();
            Dictionary<string, string> psmFilePaths = new();
            foreach(var exp in experimentFile)
            {
                string potentialPsmFilePath = Path.Combine(pathToResultsFolder, exp.SampleName, "psm.tsv");
                if(File.Exists(potentialPsmFilePath))
                {
                    psmFilePaths.Add(exp.SampleName, potentialPsmFilePath);
                }
            }

            if(psmFilePaths.Count == 0)
            {
                throw new FileNotFoundException("Could not find any PSM files in the individual file folders");
            }


            return new MsFraggerMasterFile(experimentFilePath, peptideFilePath, psmFilePaths);
        }

        public MsFraggerMasterFile(string experimentFilePath, string peptideFilePath, Dictionary<string, string> individualPsmFilePaths)
        {
            ExperimentFile = new MsFraggerExperimentFile(experimentFilePath);
            ExperimentFile.LoadResults();
            if(peptideFilePath != null)
                PeptideFile = new MsFraggerPeptideFile(peptideFilePath);


            PsmFiles = new Dictionary<MsFraggerExperiment, MsFraggerPsmFile>();
            foreach(var exp in ExperimentFile)
            {
                string psmFilePath = individualPsmFilePaths[exp.SampleName];
                PsmFiles.Add(exp, new MsFraggerPsmFile(psmFilePath));
            }
        }

        public IEnumerable<MsFraggerPsm> LoadAllPsms()
        {
            foreach(var kvp in PsmFiles)
            {
                if(kvp.Value.Count() == 0)
                    kvp.Value.LoadResults();
                foreach (var psm in kvp.Value)
                {
                    psm.SpectrumFilePath = kvp.Key.FullFilePathWithExtension;
                    yield return psm;
                }
            }
        }

    }
}
