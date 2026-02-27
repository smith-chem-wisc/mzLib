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

    public GlycoPsmFromTsv(string line, char[] split, Dictionary<string, int> parsedHeader, SpectrumMatchParsingParameters? parsingParams = null) : base(line, split, parsedHeader, parsingParams ??= new())
    {
        var spl = line.Split(split).Select(p => p.Trim('\"')).ToArray();

        GlycanMass = GetOptionalValue<double>(SpectrumMatchFromTsvHeader.GlycanMass, parsedHeader, spl);
        GlycanComposition = GetOptionalValue(SpectrumMatchFromTsvHeader.GlycanComposition, parsedHeader, spl);
        GlycanStructure = GetOptionalValue(SpectrumMatchFromTsvHeader.GlycanStructure, parsedHeader, spl);
        var localizationLevel = GetOptionalValue(SpectrumMatchFromTsvHeader.GlycanLocalizationLevel, parsedHeader, spl);
        if (localizationLevel != null)
        {
            if (localizationLevel.Equals("NA"))
                GlycanLocalizationLevel = null;
            else
                GlycanLocalizationLevel = (LocalizationLevel)Enum.Parse(typeof(LocalizationLevel), localizationLevel);
        }
        LocalizedGlycanInPeptide = GetOptionalValue(SpectrumMatchFromTsvHeader.LocalizedGlycanInPeptide, parsedHeader, spl);
        LocalizedGlycanInProtein = GetOptionalValue(SpectrumMatchFromTsvHeader.LocalizedGlycanInProtein, parsedHeader, spl);
        R138144 = GetOptionalValue<double>(SpectrumMatchFromTsvHeader.R138144, parsedHeader, spl, 0).GetValueOrDefault(0);
        LocalizedScores = GetOptionalValue<double>(SpectrumMatchFromTsvHeader.LocalizedScores, parsedHeader, spl, 0);
        YionScore = GetOptionalValue<double>(SpectrumMatchFromTsvHeader.YionScore, parsedHeader, spl);
        DiagonosticIonScore = GetOptionalValue<double>(SpectrumMatchFromTsvHeader.DiagonosticIonScore, parsedHeader, spl);
        TotalGlycanSites = GetOptionalValue(SpectrumMatchFromTsvHeader.TotalGlycanSite, parsedHeader, spl);
        FlankingResidues = GetOptionalValue(SpectrumMatchFromTsvHeader.FlankingResidues, parsedHeader, spl);
        AllPotentialGlycanLocalization = GetOptionalValue(SpectrumMatchFromTsvHeader.AllPotentialGlycanLocalization, parsedHeader, spl);
        AllSiteSpecificLocalizationProbability = GetOptionalValue(SpectrumMatchFromTsvHeader.AllSiteSpecificLocalizationProbability, parsedHeader, spl);
        NGlycanMotifCheck = GetOptionalValue(SpectrumMatchFromTsvHeader.NGlycanMotifCheck, parsedHeader, spl);
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