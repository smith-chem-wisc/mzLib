using System.Globalization;

namespace Readers;

public class GlycoPsmFromTsv : PsmFromTsv
{
    public string GlycanStructure { get; set; }
    public double? GlycanMass { get; set; }
    public string GlycanComposition { get; set; }
    public LocalizationLevel? GlycanLocalizationLevel { get; set; }
    public string LocalizedGlycanInPeptide { get; set; }
    public string LocalizedGlycanInProtein { get; set; }
    public double R138144 { get; set; }
    public double? LocalizedScores { get; set; }
    public double? YionScore { get; set; }
    public double? DiagonosticIonScore { get; set; }
    public string TotalGlycanSites { get; set; }

    public string FlankingResidues { get; set; }
    public string AllPotentialGlycanLocalization { get; set; }
    public string NGlycanMotifCheck { get; set; }
    public string AllSiteSpecificLocalizationProbability { get; set; }

    public GlycoPsmFromTsv(string line, char[] split, Dictionary<string, int> parsedHeader) : base(line, split, parsedHeader)
    {
        var spl = line.Split(split).Select(p => p.Trim('\"')).ToArray();

        GlycanMass = (parsedHeader[SpectrumMatchFromTsvHeader.GlycanMass] < 0) ? null : (double?)double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.GlycanMass]], CultureInfo.InvariantCulture);
        GlycanComposition = (parsedHeader[SpectrumMatchFromTsvHeader.GlycanComposition] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.GlycanComposition]];
        GlycanStructure = (parsedHeader[SpectrumMatchFromTsvHeader.GlycanStructure] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.GlycanStructure]];
        var localizationLevel = (parsedHeader[SpectrumMatchFromTsvHeader.GlycanLocalizationLevel] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.GlycanLocalizationLevel]];
        if (localizationLevel != null)
        {
            if (localizationLevel.Equals("NA"))
                GlycanLocalizationLevel = null;
            else
                GlycanLocalizationLevel = (LocalizationLevel)Enum.Parse(typeof(LocalizationLevel), localizationLevel);
        }
        LocalizedGlycanInPeptide = (parsedHeader[SpectrumMatchFromTsvHeader.LocalizedGlycanInPeptide] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.LocalizedGlycanInPeptide]];
        LocalizedGlycanInProtein = (parsedHeader[SpectrumMatchFromTsvHeader.LocalizedGlycanInProtein] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.LocalizedGlycanInProtein]];
        R138144 = (parsedHeader[SpectrumMatchFromTsvHeader.R138144] < 0) ? 0 : double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.R138144]], CultureInfo.InvariantCulture);
        LocalizedScores = (parsedHeader[SpectrumMatchFromTsvHeader.LocalizedScores] < 0) ? 0 : double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.LocalizedScores]], CultureInfo.InvariantCulture);
        YionScore = (parsedHeader[SpectrumMatchFromTsvHeader.YionScore] < 0) ? null : (double?)double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.YionScore]], CultureInfo.InvariantCulture);
        DiagonosticIonScore = (parsedHeader[SpectrumMatchFromTsvHeader.DiagonosticIonScore] < 0) ? null : (double?)double.Parse(spl[parsedHeader[SpectrumMatchFromTsvHeader.DiagonosticIonScore]], CultureInfo.InvariantCulture);
        TotalGlycanSites = (parsedHeader[SpectrumMatchFromTsvHeader.TotalGlycanSite] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.TotalGlycanSite]];
        FlankingResidues = (parsedHeader[SpectrumMatchFromTsvHeader.FlankingResidues] < 0) ? null
            : spl[parsedHeader[SpectrumMatchFromTsvHeader.FlankingResidues]];
        AllPotentialGlycanLocalization = (parsedHeader[SpectrumMatchFromTsvHeader.AllPotentialGlycanLocalization] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.AllPotentialGlycanLocalization]];
        AllSiteSpecificLocalizationProbability = (parsedHeader[SpectrumMatchFromTsvHeader.AllSiteSpecificLocalizationProbability] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.AllSiteSpecificLocalizationProbability]];
        NGlycanMotifCheck = (parsedHeader[SpectrumMatchFromTsvHeader.NGlycanMotifCheck] < 0) ? null : spl[parsedHeader[SpectrumMatchFromTsvHeader.NGlycanMotifCheck]];
    }

    public GlycoPsmFromTsv(GlycoPsmFromTsv psm, string fullSequence, int index = 0, string baseSequence = "") : base(psm, fullSequence, index, baseSequence)
    {
        GlycanStructure = psm.GlycanStructure;
        GlycanMass = psm.GlycanMass;
        GlycanComposition = psm.GlycanComposition;
        GlycanLocalizationLevel = psm.GlycanLocalizationLevel;
        LocalizedGlycanInPeptide = psm.LocalizedGlycanInPeptide;
        LocalizedGlycanInProtein = psm.LocalizedGlycanInProtein;
        R138144 = psm.R138144;
        LocalizedScores = psm.LocalizedScores;
        YionScore = psm.YionScore;
        DiagonosticIonScore = psm.DiagonosticIonScore;
        TotalGlycanSites = psm.TotalGlycanSites;
    }

    public static List<Tuple<int, string, double>> ReadLocalizedGlycan(string localizedGlycan)
    {
        List<Tuple<int, string, double>> tuples = new List<Tuple<int, string, double>>();
        if (localizedGlycan == null)
        {
            return tuples;
        }
        var lgs = localizedGlycan.Split(new string[] { "[", "]" }, StringSplitOptions.RemoveEmptyEntries);
        foreach (var lg in lgs)
        {
            var g = lg.Split(',', StringSplitOptions.RemoveEmptyEntries);

            Tuple<int, string, double> tuple = new Tuple<int, string, double>(int.Parse(g[0], CultureInfo.InvariantCulture), g[1], double.Parse(g[2], CultureInfo.InvariantCulture));
            tuples.Add(tuple);
        }

        return tuples;
    }   
}