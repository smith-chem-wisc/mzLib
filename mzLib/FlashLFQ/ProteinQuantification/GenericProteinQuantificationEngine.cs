using Easy.Common.Extensions;
using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace FlashLFQ
{
    /// <summary>
    /// Defines the strategy interface for protein quantification algorithms
    /// </summary>
    public interface IProteinQuantificationStrategy<TPeptide, TProteinGroup, TSpectraFile>
        where TPeptide : IQuantifiablePeptide
        where TProteinGroup : IQuantifiableProteinGroup
        where TSpectraFile : IQuantifiableSpectraFile
    {
        /// <summary>
        /// Quantifies proteins based on their constituent peptides
        /// </summary>
        void QuantifyProteins(
            Dictionary<TProteinGroup, List<TPeptide>> proteinGroupToPeptides,
            List<TSpectraFile> spectraFiles,
            bool useSharedPeptides);
    }

    /// <summary>
    /// Generic protein quantification engine that uses a strategy pattern
    /// to support different quantification algorithms (Top3, MedianPolish, etc.)
    /// </summary>
    public class GenericProteinQuantificationEngine<TPeptide, TProteinGroup, TSpectraFile>
        where TPeptide : IQuantifiablePeptide
        where TProteinGroup : IQuantifiableProteinGroup
        where TSpectraFile : IQuantifiableSpectraFile
    {
        private readonly IProteinQuantificationStrategy<TPeptide, TProteinGroup, TSpectraFile> _strategy;
        private readonly Dictionary<string, TPeptide> _peptideModifiedSequences;
        private readonly Dictionary<string, TProteinGroup> _proteinGroups;
        private readonly List<TSpectraFile> _spectraFiles;

        public GenericProteinQuantificationEngine(
            IProteinQuantificationStrategy<TPeptide, TProteinGroup, TSpectraFile> strategy,
            Dictionary<string, TPeptide> peptideModifiedSequences,
            Dictionary<string, TProteinGroup> proteinGroups,
            List<TSpectraFile> spectraFiles)
        {
            _strategy = strategy ?? throw new ArgumentNullException(nameof(strategy));
            _peptideModifiedSequences = peptideModifiedSequences ?? throw new ArgumentNullException(nameof(peptideModifiedSequences));
            _proteinGroups = proteinGroups ?? throw new ArgumentNullException(nameof(proteinGroups));
            _spectraFiles = spectraFiles ?? throw new ArgumentNullException(nameof(spectraFiles));
        }

        /// <summary>
        /// Runs the protein quantification using the configured strategy
        /// </summary>
        public void Run(bool useSharedPeptides = false)
        {
            // Reset protein intensities to 0
            foreach (var proteinGroup in _proteinGroups)
            {
                foreach (TSpectraFile file in _spectraFiles)
                {
                    proteinGroup.Value.SetIntensity(file, 0);
                }
            }

            // Associate peptides with proteins
            Dictionary<TProteinGroup, List<TPeptide>> proteinGroupToPeptides = 
                new Dictionary<TProteinGroup, List<TPeptide>>();

            List<TPeptide> peptides = _peptideModifiedSequences.Values
                .Where(p => p.UnambiguousPeptideQuant())
                .ToList();

            foreach (TPeptide peptide in peptides)
            {
                if (!peptide.UseForProteinQuant)
                {
                    continue;
                }

                // Handle shared peptides based on useSharedPeptides parameter
                bool isSharedPeptide = peptide.ProteinGroups.Count() > 1;
                if (!useSharedPeptides && isSharedPeptide)
                {
                    continue;
                }

                foreach (var pg in peptide.ProteinGroups)
                {
                    TProteinGroup proteinGroup = (TProteinGroup)pg;
                    
                    if (proteinGroupToPeptides.TryGetValue(proteinGroup, out var peptidesForThisProtein))
                    {
                        peptidesForThisProtein.Add(peptide);
                    }
                    else
                    {
                        proteinGroupToPeptides.Add(proteinGroup, new List<TPeptide> { peptide });
                    }
                }
            }

            // Execute the quantification strategy
            _strategy.QuantifyProteins(proteinGroupToPeptides, _spectraFiles, useSharedPeptides);
        }
    }
}
