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
        /// the rule MetaMorpheus will adopt for its own quantification. One confidence value decides, the first tier
        /// the file supports (<see cref="QuantifiedPsmRule.ChooseTier"/>):
        /// <list type="bullet">
        /// <item>MetaMorpheus results (<see cref="SpectrumMatchFromTsv"/>, <see cref="LightWeightSpectralMatch"/>): the PEP
        /// q-value when the file's PEP was trained, otherwise the notch q-value when the file has one, otherwise the
        /// q-value;</item>
        /// <item>DIA-NN (<see cref="DiaNnPrecursor"/>): Global.Q.Value. DIA-NN reports no PEP q-value (its PEP column is a
        /// raw posterior error probability) and no notch, so that is the only tier;</item>
        /// <item>a record type that reports no q-value (MSFragger PSMs, for one) is kept: its tool filtered it;</item>
        /// <item>an ambiguous match is never kept, whatever its type.</item>
        /// </list>
        /// The tier is decided once per result file, over its target records. Decoys are kept, as
        /// <see cref="MakeIdentifications"/> keeps them: match-between-runs needs them. Each identification's q-value is
        /// the tier's value (the q-value where a match has no notch q-value).
        /// </summary>
        /// <exception cref="MzLibUtil.MzLibException">Neither PEP nor the notch is usable and no record carries a q-value
        /// (the file has no q-value column), so every match would be dropped.</exception>
        public static List<Identification> MakeQuantifiedIdentifications(this IQuantifiableResultFile quantifiable,
            List<SpectraFileInfo> spectraFiles, double threshold = QuantifiedPsmRule.DefaultThreshold) =>
            quantifiable.MakeQuantifiedIdentifications(spectraFiles, out _, threshold);

        /// <summary>
        /// <see cref="MakeQuantifiedIdentifications(IQuantifiableResultFile, List{SpectraFileInfo}, double)"/>, also giving
        /// the <paramref name="tier"/> chosen for this file, so a caller can report which value filtered it.
        /// </summary>
        public static List<Identification> MakeQuantifiedIdentifications(this IQuantifiableResultFile quantifiable,
            List<SpectraFileInfo> spectraFiles, out QuantifiedPsmTier tier, double threshold = QuantifiedPsmRule.DefaultThreshold)
        {
            List<IQuantifiableRecord> records = quantifiable.GetQuantifiableResults().ToList();
            var targetConfidences = records.Where(record => !record.IsDecoy)
                .Select(ConfidenceOf).Where(confidence => confidence is not null).Select(confidence => confidence.Value).ToList();
            tier = QuantifiedPsmRule.ChooseTier(
                targetConfidences.Where(c => c.PepQValue is not null).Select(c => c.PepQValue.Value),
                targetConfidences.Where(c => c.NotchQValue is not null).Select(c => c.NotchQValue.Value),
                targetConfidences.Where(c => c.Pep is not null).Select(c => c.Pep.Value));
            if (tier == QuantifiedPsmTier.QValue && targetConfidences.Count > 0
                && targetConfidences.All(confidence => double.IsNaN(confidence.QValue)))
            {
                throw new MzLibUtil.MzLibException("No match in the result file has a q-value (is its q-value column missing?), " +
                    "and it has no usable PEP q-value or notch q-value, so no match could be quantified.");
            }

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

                if (ConfidenceOf(record) is var (qValue, notchQValue, pepQValue, _)
                    && !QuantifiedPsmRule.PassesConfidence(tier, qValue, notchQValue, pepQValue, threshold))
                {
                    continue;
                }

                identifications.Add(MakeIdentification(record, spectraFile, allProteinGroups, tier));
            }

            return identifications;
        }

        /// <summary>
        /// The confidence values <see cref="QuantifiedPsmRule"/> reads from a record, or null for a record type that
        /// reports none. PEP is read only to tell whether PEP was trained.
        /// </summary>
        private static (double QValue, double? NotchQValue, double? PepQValue, double? Pep)? ConfidenceOf(IQuantifiableRecord record) => record switch
        {
            SpectrumMatchFromTsv psm => (psm.QValue, psm.QValueNotch, psm.PEP_QValue, psm.PEP),
            LightWeightSpectralMatch light => (light.QValue, light.QValueNotch, light.PepQValue, light.Pep),
            DiaNnPrecursor diaNn => (diaNn.GlobalQValue, null, null, null),
            _ => null,
        };

        private static Identification MakeIdentification(IQuantifiableRecord record, SpectraFileInfo spectraFile,
            Dictionary<string, ProteinGroup> allProteinGroups, bool usePepQValue) =>
            MakeIdentification(record, spectraFile, allProteinGroups, usePepQValue ? QuantifiedPsmTier.PepQValue : QuantifiedPsmTier.QValue);

        /// <summary>
        /// One record as a FlashLFQ identification, its protein groups shared through <paramref name="allProteinGroups"/>.
        /// Its q-value is the one <paramref name="tier"/> names.
        /// </summary>
        private static Identification MakeIdentification(IQuantifiableRecord record, SpectraFileInfo spectraFile,
            Dictionary<string, ProteinGroup> allProteinGroups, QuantifiedPsmTier tier)
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
            double? notchQValue = null;
            double score = 0;
            // Populate optional fields, which only some result types report
            if( record is SpectrumMatchFromTsv psmFromTsv)
            {
                qValue = psmFromTsv.QValue;
                pepQValue = psmFromTsv.PEP_QValue;
                notchQValue = psmFromTsv.QValueNotch;
                score = psmFromTsv.Score;
            }
            else if (record is LightWeightSpectralMatch light)
            {
                qValue = light.QValue;
                pepQValue = light.PepQValue;
                notchQValue = light.QValueNotch;
                score = light.Score;
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
                qValue: tier switch
                {
                    QuantifiedPsmTier.PepQValue => pepQValue,
                    QuantifiedPsmTier.QValueNotch => notchQValue ?? qValue,
                    _ => qValue,
                });
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