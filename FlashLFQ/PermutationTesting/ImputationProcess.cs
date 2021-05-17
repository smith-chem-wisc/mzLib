using MathNet.Numerics.Distributions;
using System;
using System.Collections.Generic;
using System.Linq;

namespace FlashLFQ
{
    /// <summary>
    /// Class used for imputing missing intensity values for each protein
    /// </summary>
    public class ImputationProcess
    {
        /// <summary>
        /// Imputes missing intensity values
        /// </summary>
        public List<ProteinGroup> RunImputationProcess(List<ProteinGroup> proteins, List<SpectraFileInfo> sampleFileNames, double meanFraction)
        {
            //foreach sample, determine the distribution of protein intensities
            //create a function to impute data for the given sample
            //impute missing values

            //create array to store all intensities for each file
            List<double>[] intensitiesForEachSample = new List<double>[sampleFileNames.Count];
            for(int i=0; i<sampleFileNames.Count; i++)
            {
                intensitiesForEachSample[i] = new List<double>();
            }

            //populate intensities for each file
            foreach(ProteinGroup protein in proteins)
            {
                for(int i=0; i<sampleFileNames.Count; i++)
                {
                   intensitiesForEachSample[i].Add(protein.GetIntensity(sampleFileNames[i]));
                }
            }

            //create sample-specific functions for imputation
            Normal[] imputationFunctions = new Normal[sampleFileNames.Count];
            for(int i=0; i<sampleFileNames.Count; i++)
            {
                List<double> intensitiesForThisSample = intensitiesForEachSample[i];
                double mean = intensitiesForThisSample.Average();
                double stddev = Math.Sqrt((intensitiesForThisSample.Sum(x => Math.Pow(x - mean, 2))) / intensitiesForThisSample.Count);
                double fractionOfMissingValues = (proteins.Count - intensitiesForThisSample.Count) * 1d / proteins.Count;
                imputationFunctions[i] = CreateImputationFunction(mean, stddev, fractionOfMissingValues, meanFraction);
            }

            //impute values for each protein
            List<ProteinGroup> proteinsWithImputedValues = new List<ProteinGroup>(proteins.Count);
            foreach(ProteinGroup protein in proteins)
            {
                //create new protein to store the imputed values
                ProteinGroup imputedProtein = new ProteinGroup(protein.ProteinGroupName, protein.GeneName, protein.Organism);
                for(int i=0; i<sampleFileNames.Count; i++)
                {
                    SpectraFileInfo file = sampleFileNames[i];
                    double intensity = protein.GetIntensity(file);
                    intensity = intensity > 0 ? intensity : imputationFunctions[i].Sample(); //impute if needed
                    imputedProtein.SetIntensity(file, intensity);
                }
                proteinsWithImputedValues.Add(imputedProtein);
            }

            return proteinsWithImputedValues;
        }

        /// <summary>
        /// Imputes missing intensity value for each protein
        /// </summary>
        public Normal CreateImputationFunction(double mean, double stddev, double fractionOfMissingValues, double meanFraction)
        {
            double imputedProbability = Math.Min(fractionOfMissingValues / (1 - fractionOfMissingValues), 1); //TODO: Explanation needed
            double standardDeviationFraction = Math.Max(2 * fractionOfMissingValues, 0.3); //TODO: Explanation needed
            double stdDevFraction = 0.6 * (1 - (fractionOfMissingValues * fractionOfMissingValues));//TODO: Explanation needed
            Normal probabilityDist = new Normal(mean, standardDeviationFraction);
            double probabilitySetPoint = probabilityDist.Density(mean + stdDevFraction * standardDeviationFraction);
            double yCoordinate = imputedProbability * probabilitySetPoint;
            double deltaX = standardDeviationFraction * stdDevFraction;
            Normal xCoord = new Normal(mean, stddev);
            double deltaMu = xCoord.InverseCumulativeDistribution(yCoordinate);
            double meanDownshift = (deltaMu - deltaX * meanFraction);

            return new Normal(meanDownshift, standardDeviationFraction, new Random(2));
        }
    }
}