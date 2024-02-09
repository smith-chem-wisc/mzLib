using System.Globalization;
using Omics.SpectrumMatch;

namespace Transcriptomics
{
    public class OsmFromTsv : SpectrumMatchFromTsv
    {
        public OsmFromTsv(string line, char[] split, Dictionary<string, int> parsedHeader)
        {
            var spl = line.Split(split).Select(p => p.Trim('\"')).ToArray();

            //Required properties
            FileNameWithoutExtension = spl[parsedHeader[SpectrumMatchFromTsvHeader.FileName]].Trim();

            // remove file format, e.g., .raw, .mzML, .mgf
            // this is more robust but slower than Path.GetFileNameWithoutExtension
            if (FileNameWithoutExtension.Contains('.'))
            {
                foreach (var knownSpectraFileExtension in SpectrumMatchFromTsvHeader.AcceptedSpectraFormats)
                {
                    FileNameWithoutExtension = Path.GetFileName(FileNameWithoutExtension.Replace(knownSpectraFileExtension, string.Empty, StringComparison.InvariantCultureIgnoreCase));
                }
            }

            Ms2ScanNumber = int.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.Ms2ScanNumber]]);

            // this will probably not be known in an .mgf data file
            if (int.TryParse(spl[parsedHeader[SpectrumMatchFromTsvHeader.PrecursorScanNum]].Trim(), out int result))
            {
                PrecursorScanNum = result;
            }
            else
            {
                PrecursorScanNum = 0;
            }

            PrecursorCharge = (int)double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.PrecursorCharge]].Trim(), CultureInfo.InvariantCulture);
            PrecursorMz = double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.PrecursorMz]].Trim(), CultureInfo.InvariantCulture);
            PrecursorMass = double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.PrecursorMass]].Trim(), CultureInfo.InvariantCulture);
            BaseSeq = RemoveParentheses(spl[parsedHeader[SpectrumMatchFromTsvHeader.BaseSequence]].Trim());
            FullSequence = spl[parsedHeader[SpectrumMatchFromTsvHeader.FullSequence]];
            MonoisotopicMass = spl[parsedHeader[SpectrumMatchFromTsvHeader.MonoisotopicMass]].Trim();
            Score = double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.Score]].Trim(), CultureInfo.InvariantCulture);
            DecoyContamTarget = spl[parsedHeader[SpectrumMatchFromTsvHeader.DecoyContaminantTarget]].Trim();

            //we are reading in all primary and child ions here only to delete the child scans later. This should be done better.
            MatchedIons = (spl[parsedHeader[SpectrumMatchFromTsvHeader.MatchedIonMzRatios]].StartsWith("{")) ?
                ReadChildScanMatchedIons(spl[parsedHeader[SpectrumMatchFromTsvHeader.MatchedIonMzRatios]].Trim(), spl[parsedHeader[SpectrumMatchFromTsvHeader.MatchedIonIntensities]].Trim(), BaseSeq).First().Value :
                ReadFragmentIonsFromString(spl[parsedHeader[SpectrumMatchFromTsvHeader.MatchedIonMzRatios]].Trim(), spl[parsedHeader[SpectrumMatchFromTsvHeader.MatchedIonIntensities]].Trim(), BaseSeq, spl[parsedHeader[SpectrumMatchFromTsvHeader.MatchedIonMassDiffDa]].Trim());

            AmbiguityLevel = (parsedHeader[SpectrumMatchFromTsvHeader.AmbiguityLevel] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.AmbiguityLevel]].Trim();

            QValue = double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.QValue]].Trim(), CultureInfo.InvariantCulture);
            //For general psms
            TotalIonCurrent = (parsedHeader[SpectrumMatchFromTsvHeader.TotalIonCurrent] < 0) ? null : (double?)double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.TotalIonCurrent]].Trim(), CultureInfo.InvariantCulture);
            DeltaScore = (parsedHeader[SpectrumMatchFromTsvHeader.DeltaScore] < 0) ? null : (double?)double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.DeltaScore]].Trim(), CultureInfo.InvariantCulture);
            Notch = (parsedHeader[SpectrumMatchFromTsvHeader.Notch] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.Notch]].Trim();
            EssentialSeq = (parsedHeader[SpectrumMatchFromTsvHeader.EssentialSequence] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.EssentialSequence]].Trim();
            MissedCleavage = (parsedHeader[SpectrumMatchFromTsvHeader.MissedCleavages] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.MissedCleavages]].Trim();
            MassDiffDa = (parsedHeader[SpectrumMatchFromTsvHeader.MassDiffDa] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.MassDiffDa]].Trim();
            MassDiffPpm = (parsedHeader[SpectrumMatchFromTsvHeader.MassDiffPpm] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.MassDiffPpm]].Trim();
            Accession = (parsedHeader[SpectrumMatchFromTsvHeader.Accession] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.Accession]].Trim();
            Name = (parsedHeader[SpectrumMatchFromTsvHeader.Name] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.Name]].Trim();
            GeneName = (parsedHeader[SpectrumMatchFromTsvHeader.GeneName] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.GeneName]].Trim();
            OrganismName = (parsedHeader[SpectrumMatchFromTsvHeader.OrganismName] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.OrganismName]].Trim();
            IntersectingSequenceVariations = (parsedHeader[SpectrumMatchFromTsvHeader.IntersectingSequenceVariations] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.IntersectingSequenceVariations]].Trim();
            IdentifiedSequenceVariations = (parsedHeader[SpectrumMatchFromTsvHeader.IdentifiedSequenceVariations] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.IdentifiedSequenceVariations]].Trim();
            SpliceSites = (parsedHeader[SpectrumMatchFromTsvHeader.SpliceSites] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.SpliceSites]].Trim();
            Description = (parsedHeader[SpectrumMatchFromTsvHeader.Description] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.Description]].Trim();
            StartAndEndResiduesInParentSequence = (parsedHeader[SpectrumMatchFromTsvHeader.StartAndEndResiduesInFullSequence] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.StartAndEndResiduesInFullSequence]].Trim();
            PreviousResidue = (parsedHeader[SpectrumMatchFromTsvHeader.PreviousResidue] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.PreviousResidue]].Trim();
            NextResidue = (parsedHeader[SpectrumMatchFromTsvHeader.NextResidue] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.NextResidue]].Trim();
            QValueNotch = (parsedHeader[SpectrumMatchFromTsvHeader.QValueNotch] < 0) ? null : (double?)double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.QValueNotch]].Trim(), CultureInfo.InvariantCulture);
            RetentionTime = (parsedHeader[SpectrumMatchFromTsvHeader.Ms2ScanRetentionTime] < 0) ? null : (double?)double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.Ms2ScanRetentionTime]].Trim(), CultureInfo.InvariantCulture);
            PEP = double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.PEP]].Trim(), CultureInfo.InvariantCulture);
            PEP_QValue = double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.PEP_QValue]].Trim(), CultureInfo.InvariantCulture);
            VariantCrossingIons = FindVariantCrossingIons();
            SpectralAngle = (parsedHeader[SpectrumMatchFromTsvHeader.SpectralAngle] < 0) ? null : (double?)double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.SpectralAngle]].Trim(), CultureInfo.InvariantCulture);

        }
    }
}
