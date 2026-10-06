using System;
using System.Collections.Generic;
using System.Linq;
using Readers;
using MassSpectrometry;
using Quantification;

namespace FlashLFQ
{
    public static class MzLibExtensions
    {
        /// <summary>
        /// Makes a list of identification objects usable by FlashLFQ from an IQuantifiableResultFile
        /// </summary>
        public static List<Identification> MakeIdentifications(this IQuantifiableResultFile quantifiable, List<SpectraFileInfo> spectraFiles, bool usePepQValue = false)
        {
            IEnumerable<IQuantifiableRecord> quantifiableRecords = quantifiable.GetQuantifiableResults();
            List<Identification> identifications = new List<Identification>();
            Dictionary<string, ProteinGroup> allProteinGroups = new Dictionary<string, ProteinGroup>();
            Dictionary<string, SpectraFileInfo> allSpectraFiles = MakeSpectraFileDict(quantifiable, spectraFiles);

            foreach (var record in quantifiableRecords)
            {
                // Get the spectra file info from the dictionary using the file name. A result file can
                // reference more spectra files than are supplied for quantification (e.g. an .osmtsv or
                // .psmtsv covering runs the user did not load). Skip identifications for any file that was
                // not supplied rather than aborting the whole read - this matches the legacy PsmReader path (in the FlashLFQ Standalone),
                // which returns null for PSMs whose spectrum file has no data input.
                if (!allSpectraFiles.TryGetValue(record.FileName, out var spectraFile) || spectraFile is null)
                {
                    continue;
                }

                identifications.Add(MakeIdentification(record, spectraFile, allProteinGroups, usePepQValue));
            }

            return identifications;
        }

        /// <summary>
        /// Like <see cref="MakeIdentifications"/>, but keeps only the matches <see cref="QuantifiedPsmRule"/> quantifies,
        /// the same rule MetaMorpheus applies to its own quantification:
        /// <list type="bullet">
        /// <item>MetaMorpheus results (<see cref="SpectrumMatchFromTsv"/>, <see cref="LightWeightSpectralMatch"/>): the PEP
        /// q-value when the file's PEP was trained, otherwise the q-value and the notch q-value (where the file has one);</item>
        /// <item>DIA-NN (<see cref="DiaNnPrecursor"/>): Global.Q.Value. DIA-NN reports no PEP q-value (its PEP column is a
        /// raw posterior error probability) and no notch, so that is the only tier;</item>
        /// <item>a record type that reports no q-value (MSFragger PSMs, for one) is kept: its tool filtered it;</item>
        /// <item>an ambiguous match is never kept, whatever its type.</item>
        /// </list>
        /// Whether PEP was trained is decided once, over every record in the file. Decoys are kept, as
        /// <see cref="MakeIdentifications"/> keeps them: match-between-runs needs them.
        /// </summary>
        public static List<Identification> MakeQuantifiedIdentifications(this IQuantifiableResultFile quantifiable,
            List<SpectraFileInfo> spectraFiles, double threshold = QuantifiedPsmRule.DefaultThreshold)
        {
            List<IQuantifiableRecord> records = quantifiable.GetQuantifiableResults().ToList();
            bool usePepQValue = QuantifiedPsmRule.PepQValueIsUsable(records
                .Select(ConfidenceOf)
                .Where(confidence => confidence?.PepQValue is not null)
                .Select(confidence => confidence.Value.PepQValue.Value));

            List<Identification> identifications = new List<Identification>();
            Dictionary<string, ProteinGroup> allProteinGroups = new Dictionary<string, ProteinGroup>();
            Dictionary<string, SpectraFileInfo> allSpectraFiles = MakeSpectraFileDict(quantifiable, spectraFiles);

            foreach (var record in records)
            {
                if (!allSpectraFiles.TryGetValue(record.FileName, out var spectraFile) || spectraFile is null)
                {
                    continue;
                }

                if (QuantifiedPsmRule.IsAmbiguous(record.BaseSequence, record.FullSequence))
                {
                    continue;
                }

                if (ConfidenceOf(record) is var (qValue, notchQValue, pepQValue)
                    && !QuantifiedPsmRule.PassesConfidence(qValue, notchQValue, pepQValue, usePepQValue, threshold))
                {
                    continue;
                }

                identifications.Add(MakeIdentification(record, spectraFile, allProteinGroups, usePepQValue));
            }

            return identifications;
        }

        /// <summary>
        /// The confidence values <see cref="QuantifiedPsmRule"/> reads from a record, or null for a record type that
        /// reports none.
        /// </summary>
        private static (double QValue, double? NotchQValue, double? PepQValue)? ConfidenceOf(IQuantifiableRecord record) => record switch
        {
            SpectrumMatchFromTsv psm => (psm.QValue, psm.QValueNotch, psm.PEP_QValue),
            LightWeightSpectralMatch light => (light.QValue, null, light.PepQValue),
            DiaNnPrecursor diaNn => (diaNn.GlobalQValue, null, null),
            _ => null,
        };

        /// <summary>
        /// One record as a FlashLFQ identification, its protein groups shared through <paramref name="allProteinGroups"/>.
        /// </summary>
        private static Identification MakeIdentification(IQuantifiableRecord record, SpectraFileInfo spectraFile,
            Dictionary<string, ProteinGroup> allProteinGroups, bool usePepQValue)
        {
            string baseSequence = record.BaseSequence;
            string modifiedSequence = record.FullSequence;
            double ms2RetentionTimeInMinutes = record.RetentionTime;
            double monoisotopicMass = record.MonoisotopicMass;
            int precursorChargeState = record.ChargeState;

            List<ProteinGroup> proteinGroups = new();
            foreach (var info in record.ProteinGroupInfos)
            {
                if (allProteinGroups.TryGetValue(info.proteinAccessions, out var proteinGroup))
                {
                    proteinGroups.Add(proteinGroup);
                }
                else
                {
                    allProteinGroups.Add(info.proteinAccessions, new ProteinGroup(info.proteinAccessions, info.geneName, info.organism));
                    proteinGroups.Add(allProteinGroups[info.proteinAccessions]);
                }
            }

            double qValue = 0;
            double pepQValue = 0;
            double score = 0;
            // Populate optional fields, which only some result types report
            if( record is SpectrumMatchFromTsv psmFromTsv)
            {
                qValue = psmFromTsv.QValue;
                pepQValue = psmFromTsv.PEP_QValue;
                score = psmFromTsv.Score;
            }
            else if (record is DiaNnPrecursor diaNnPrecursor)
            {
                // Global.Q.Value, not Q.Value: DIA-NN's Q.Value is scoped to a single run and is
                // already filtered below 1%, so it would let every precursor through the
                // experiment-wide gates below. Global.Q.Value is the run-spanning figure, which
                // is what DonorQValueThreshold is comparing against when it picks MBR donors.
                qValue = diaNnPrecursor.GlobalQValue;

                // DIA-NN's PEP column is the raw per-precursor posterior error probability, i.e.
                // the local probability that this one identification is wrong. It is stored as-is,
                // unlike MetaMorpheus's PEP_QValue, which is a q-value derived from PEP (a monotone
                // cumulative FDR). They are different quantities on different scales, so gating on
                // this with usePepQValue is not equivalent to gating on a PEP-derived q-value.
                pepQValue = diaNnPrecursor.PosteriorErrorProbability;

                // CScore is DIA-NN's classifier score, higher being better, which matches how
                // DonorCriterion.Score reads PsmScore
                score = diaNnPrecursor.CScore ?? 0;
            }

            return new Identification(
                spectraFile, 
                baseSequence, 
                modifiedSequence, 
                monoisotopicMass, 
                ms2RetentionTimeInMinutes, 
                precursorChargeState, 
                proteinGroups, 
                useForProteinQuant: !record.IsDecoy, 
                decoy: record.IsDecoy,
                psmScore: score,
                qValue: usePepQValue ? pepQValue : qValue);
        }

        private static Dictionary<string, SpectraFileInfo> MakeSpectraFileDict(this IQuantifiableResultFile quantifiable, List<SpectraFileInfo> spectraFiles)
        {
            Dictionary<string, SpectraFileInfo> allSpectraFiles = new Dictionary<string, SpectraFileInfo>();

            // 1. from list of SFIs create a list of strings that contains each full file path w/ extension
            List<string> fullFilePaths = new List<string>();
            foreach (SpectraFileInfo spectraFileInfo in spectraFiles)
            {
                fullFilePaths.Add(spectraFileInfo.FullFilePathWithExtension);
            }

            // 2. call quantifiableresultfile.filenametofilepath and get stringstring dict
            Dictionary<string, string> allFiles = quantifiable.FileNameToFilePath(fullFilePaths);

            // 3. using stringstring dict create a string spectrafileinfo dict where key is same b/w dicts and value fullfilepath is replaced spectrafileinfo obj
            foreach (var file in allFiles)
            {
                string key = file.Key;
                string filePath = file.Value;
                // FirstOrDefault matches the 1st elt from spectraFiles w/ specified filePath
                SpectraFileInfo? matchingSpectraFile = spectraFiles.FirstOrDefault(spectraFileInfo => spectraFileInfo.FullFilePathWithExtension == filePath);
                if (!allSpectraFiles.ContainsKey(key))
                {
                    allSpectraFiles[key] = matchingSpectraFile;
                }
            }

            return allSpectraFiles;
        }
    }
}