using MassSpectrometry;
using Omics;
using Omics.BioPolymerGroup;
using Quantification.Interfaces;

namespace Quantification.Strategies
{
    public abstract class RollUpStrategyBase : IRollUpStrategy
    {
        public abstract string Name { get; }
        public abstract PeptideMatrix RollUpSpectralMatches(IExperimentalDesign experimentalDesign, List<ISpectralMatch> spectralMatches, List<IBioPolymerWithSetMods> peptides);
        public abstract ProteinMatrix RollUpPeptides(PeptideMatrix peptides, List<IBioPolymerGroup> proteins);


        // Can get more generic
        // What we need is a function to group lower level entities (spectral matches or peptides) by higher level entities (peptides or proteins)
        // Then, a function that takes in a group of lower level entities and outputs a rolled up value for the higher level entity
        // Basically, will take in a pxn matrix and output a qxn matrix where q <= p


        /// <summary>
        /// Orders first by file name (alphabetically), then by sample info for each file
        /// </summary>
        /// <param name="experimentalDesign"></param>
        /// <param name="spectralMatches"></param>
        /// <param name="filePathToArrayPositionDict"> Dictionary mapping file paths to their corresponding zero-indexed array positions in the ordered sample list </param>
        /// <returns></returns>
        public List<ISampleInfo> GetOrderedSampleInfos(IExperimentalDesign experimentalDesign, List<ISpectralMatch> spectralMatches, 
            out Dictionary<string, List<int>> filePathToArrayPositionDict)
        {
            // Start by ordering the file paths alphabetically        
            var filePaths = spectralMatches
                .Select(sm => sm.FullFilePath)
                .Distinct()
                .OrderBy(fp => fp)
                .ToList();

            List<ISampleInfo> orderedSamples = new();
            filePathToArrayPositionDict = new Dictionary<string, List<int>>();
            int zeroIndexedArrayPosition = 0;
            foreach (var filePath in filePaths)
            {
                var fileName = Path.GetFileName(filePath);
                filePathToArrayPositionDict[filePath] = new List<int>();
                if (experimentalDesign.FileNameSampleInfoDictionary.TryGetValue(fileName, out var sampleInfos))
                {
                    orderedSamples.AddRange(sampleInfos);
                    filePathToArrayPositionDict[filePath].AddRange(Enumerable.Range(zeroIndexedArrayPosition, sampleInfos.Length));
                    zeroIndexedArrayPosition += sampleInfos.Length;
                }
                else
                {
                    throw new KeyNotFoundException($"File name '{fileName}' not found in experimental design.");
                }
            }

            return orderedSamples;
        }
    }
}
