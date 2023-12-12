﻿using MzLibUtil;
using Omics.SpectrumMatch;
using Proteomics.PSM;

namespace Readers
{
    public static class SpectrumMatchTsvReader
    {
        private static readonly char[] Split = { '\t' };

        /// <summary>
        /// Legacy method for reading PsmFromTsv files, creates a generic SpectrumMatchFromTsv object for each line
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="warnings"></param>
        /// <returns></returns>
        /// <exception cref="MzLibException"></exception>
        /// <exception cref="ArgumentOutOfRangeException"></exception>
        public static List<SpectrumMatchFromTsv> ReadTsv(string filePath, out List<string> warnings)
        {
            List<SpectrumMatchFromTsv> psms = new List<SpectrumMatchFromTsv>();
            warnings = new List<string>();

            StreamReader reader = null;
            try
            {
                reader = new StreamReader(filePath);
            }
            catch (Exception e)
            {
                throw new MzLibException("Could not read file: " + e.Message);
            }

            int lineCount = 0;

            string line;
            Dictionary<string, int> parsedHeader = null;

            var fileType = filePath.ParseFileType();
            while (reader.Peek() > 0)
            {
                lineCount++;

                line = reader.ReadLine();

                if (lineCount == 1)
                {
                    parsedHeader = ParseHeader(line);
                    continue;
                }

                try
                {
                    switch (filePath.ParseFileType())
                    {
                        case SupportedFileType.psmtsv:
                            psms.Add(new PsmFromTsv(line, Split, parsedHeader));
                            break;

                        // TODO: Create an osmtsv case when transcriptomics is merged

                        default:
                            throw new ArgumentOutOfRangeException();
                    }
                }
                catch (Exception e)
                {
                    warnings.Add("Could not read line: " + lineCount);
                }
            }

            reader.Close();

            if (lineCount - 1 != psms.Count)
            {
                warnings.Add("Warning: " + (lineCount - 1 - psms.Count) + " PSMs were not read.");
            }

            return psms;
        }

        /// <summary>
        /// Reads a psmtsv file and returns PsmFromTsv objects
        /// It is simply a cast of the ReadTsv method
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="warnings"></param>
        /// <returns></returns>
        public static List<PsmFromTsv> ReadPsmTsv(string filePath, out List<string> warnings) =>
            ReadTsv(filePath, out warnings).Cast<PsmFromTsv>().ToList();

        //TODO: Add a ReadOsmTsv method when transcriptomics is merged

        public static Dictionary<string, int> ParseHeader(string header)
        {
            var parsedHeader = new Dictionary<string, int>();
            var spl = header.Split(Split);

            parsedHeader.Add(SpectrumMatchFromTsvHeader.FullSequence, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.FullSequence));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.Ms2ScanNumber, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.Ms2ScanNumber));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.FileName, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.FileName));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.TotalIonCurrent, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.TotalIonCurrent));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.PrecursorScanNum, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.PrecursorScanNum));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.PrecursorCharge, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.PrecursorCharge));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.PrecursorMz, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.PrecursorMz));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.PrecursorMass, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.PrecursorMass));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.Score, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.Score));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.DeltaScore, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.DeltaScore));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.Notch, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.Notch));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.BaseSequence, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.BaseSequence));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.EssentialSequence, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.EssentialSequence));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.AmbiguityLevel, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.AmbiguityLevel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.MissedCleavages, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.MissedCleavages));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.MassDiffDa, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.MassDiffDa));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.MassDiffPpm, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.MassDiffPpm));

            //Handle legacy input
            if (spl.Contains(SpectrumMatchFromTsvHeader.Accession))
            {
                parsedHeader.Add(SpectrumMatchFromTsvHeader.SpectrumMatchCount, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.SpectrumMatchCount));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.MonoisotopicMass, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.MonoisotopicMass));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.Accession, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.Accession));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.Name, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.Name));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.Description, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.Description));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.StartAndEndResiduesInFullSequence, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.StartAndEndResiduesInFullSequence));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.NextResidue, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.PreviousResidue));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.PreviousResidue, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.NumExperimentalPeaks));
            }
            else
            {
                parsedHeader.Add(SpectrumMatchFromTsvHeader.SpectrumMatchCount, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.PsmCount));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.MonoisotopicMass, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.PeptideMonoMass));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.Accession, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.ProteinAccession));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.Name, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.ProteinName));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.Description, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.PeptideDescription));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.StartAndEndResiduesInFullSequence, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.StartAndEndResiduesInProtein));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.NextResidue, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.PreviousAminoAcid));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.PreviousResidue, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.NextAminoAcid));
            }

            parsedHeader.Add(SpectrumMatchFromTsvHeader.GeneName, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.GeneName));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.OrganismName, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.OrganismName));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.IntersectingSequenceVariations, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.IntersectingSequenceVariations));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.IdentifiedSequenceVariations, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.IdentifiedSequenceVariations));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.SpliceSites, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.SpliceSites));

            parsedHeader.Add(SpectrumMatchFromTsvHeader.DecoyContaminantTarget, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.DecoyContaminantTarget));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.MatchedIonMzRatios, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.MatchedIonMzRatios));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.MatchedIonIntensities, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.MatchedIonIntensities));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.MatchedIonMassDiffDa, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.MatchedIonMassDiffDa));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.SpectralAngle, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.SpectralAngle));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.QValue, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.QValue));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.QValueNotch, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.QValueNotch));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.PEP, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.PEP));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.PEP_QValue, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.PEP_QValue));

            parsedHeader.Add(SpectrumMatchFromTsvHeader.CrossTypeLabel, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.CrossTypeLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.LinkResiduesLabel, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.LinkResiduesLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.ProteinLinkSiteLabel, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.ProteinLinkSiteLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.RankLabel, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.RankLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.BetaPeptideProteinAccessionLabel, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.BetaPeptideProteinAccessionLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.BetaPeptideProteinLinkSiteLabel, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.BetaPeptideProteinLinkSiteLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.BetaPeptideBaseSequenceLabel, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.BetaPeptideBaseSequenceLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.BetaPeptideFullSequenceLabel, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.BetaPeptideFullSequenceLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.BetaPeptideTheoreticalMassLabel, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.BetaPeptideTheoreticalMassLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.BetaPeptideScoreLabel, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.BetaPeptideScoreLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.BetaPeptideRankLabel, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.BetaPeptideRankLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonsLabel, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonsLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonIntensitiesLabel, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonIntensitiesLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.XLTotalScoreLabel, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.XLTotalScoreLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.ParentIonsLabel, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.ParentIonsLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.Ms2ScanRetentionTime, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.Ms2ScanRetentionTime));


            parsedHeader.Add(SpectrumMatchFromTsvHeader.GlycanMass, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.GlycanMass));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.GlycanStructure, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.GlycanStructure));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.GlycanComposition, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.GlycanComposition));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.GlycanLocalizationLevel, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.GlycanLocalizationLevel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.LocalizedGlycan, Array.IndexOf(spl, SpectrumMatchFromTsvHeader.LocalizedGlycan));

            return parsedHeader;
        }
    }
}
