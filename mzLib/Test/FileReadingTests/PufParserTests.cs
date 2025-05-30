using System.IO;
using System.Linq;
using NUnit.Framework;
using Readers.Puf;

namespace Test.FileReadingTests;

[TestFixture]
public class PufParserTests
{
    private const string TestFile = @"DataFiles\target955.puf";
    private PufDataSet _dataSet;

    [SetUp]
    public void SetUp()
    {
        var combinedLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, TestFile);
        Assert.That(File.Exists(combinedLocation), $"Test file not found: {combinedLocation}");
        _dataSet = PufParser.Parse(combinedLocation);
    }

    [Test]
    public void CanParseDataSetVersion()
    {
        Assert.That(_dataSet, Is.Not.Null);
        Assert.That(_dataSet.Version, Is.EqualTo("1.1"));
    }

    [Test]
    public void CanParseExperiment()
    {
        Assert.That(_dataSet.Experiments, Has.Count.EqualTo(1));
        var exp = _dataSet.Experiments[0];
        Assert.That(exp.Id, Is.EqualTo("955"));
        Assert.That(exp.Comment, Does.Contain("Characterization score: 221.13"));
    }

    [Test]
    public void CanParseInstrumentData()
    {
        var exp = _dataSet.Experiments[0];
        var inst = exp.InstrumentData;
        Assert.That(inst, Is.Not.Null);
        Assert.That(inst.FragmentationMethod, Is.EqualTo("HCD"));
        Assert.That(inst.IonType, Is.EqualTo("BY"));
        Assert.That(inst.Intacts, Has.Count.EqualTo(1));
        Assert.That(inst.Fragments, Is.Not.Empty);
        Assert.That(inst.Intacts[0].MassMonoisotopic, Is.EqualTo(16153.783).Within(1e-6));
    }

    [Test]
    public void CanParseAnalysisAndSearchParameters()
    {
        var exp = _dataSet.Experiments[0];
        Assert.That(exp.Analyses, Has.Count.EqualTo(1));
        var analysis = exp.Analyses[0];
        Assert.That(analysis.Id, Is.EqualTo("1"));
        Assert.That(analysis.Type, Is.EqualTo("absolute_mass"));

        var search = analysis.SearchParameters;
        Assert.That(search, Is.Not.Null);
        Assert.That(search.Database, Is.Not.Null);
        Assert.That(search.Database.InternalName, Is.EqualTo("pseudomonas_aeruginosa_2012_11_top_down_complex"));
        Assert.That(search.PtmList, Does.Contain("Pyroglutamic acid from glutamic acid"));
        Assert.That(search.IntactTolerance, Is.EqualTo("1000"));
        Assert.That(search.IntactToleranceUnit, Is.EqualTo("Da"));
        Assert.That(search.FragmentTolerance, Is.EqualTo("15"));
        Assert.That(search.FragmentToleranceUnit, Is.EqualTo("ppm"));
        Assert.That(search.DeltaM, Is.True);
        Assert.That(search.Multiplexing, Is.False);
        Assert.That(search.ITraq, Is.False);
        Assert.That(search.Disulfide, Is.False);
        Assert.That(search.HitCriteria, Is.Not.Null);
        Assert.That(search.HitCriteria.Max, Is.EqualTo(25));
        Assert.That(search.HitCriteria.MinimumMatchesNum, Is.EqualTo(4));
    }

    [Test]
    public void CanParseResultsAndHitList()
    {
        var exp = _dataSet.Experiments[0];
        var results = exp.Analyses[0].Results;
        Assert.That(results, Is.Not.Null);
        Assert.That(results.HitLists, Has.Count.EqualTo(1));
        var hitList = results.HitLists[0];
        Assert.That(hitList.IntactId, Is.EqualTo("1"));
        Assert.That(hitList.Hits, Has.Count.EqualTo(1));
        var hit = hitList.Hits[0];
        Assert.That(hit.Id, Is.EqualTo("1"));
        Assert.That(hit.Description, Does.Contain("30S ribosomal protein S6"));
        Assert.That(hit.Signalp, Is.False);
        Assert.That(hit.Propep, Is.False);
        Assert.That(hit.SequenceLength, Is.EqualTo(139));
        Assert.That(hit.TheoreticalMass, Is.EqualTo(16154.70955469).Within(1e-8));
        Assert.That(hit.MassDifferenceDa, Is.EqualTo(-0.926554689971454).Within(1e-10));
        Assert.That(hit.MassDifferencePpm, Is.EqualTo(-57.3550819242343).Within(1e-10));
        Assert.That(hit.Score, Is.Not.Null);
        Assert.That(hit.Score.PScore, Is.EqualTo("4.99e-012"));
        Assert.That(hit.Score.Expected, Is.EqualTo("1.51e-007"));
        Assert.That(hit.Score.McluckeyScore, Is.EqualTo("0.0"));
        Assert.That(hit.MatchingSequence, Is.Not.Null);
        Assert.That(hit.MatchingSequence.Resid, Does.StartWith("MRHYEIVFLVHPDQSEQV"));
        Assert.That(hit.MatchingFragments, Is.Not.Empty);
        Assert.That(results.Legend, Is.EqualTo("Cysteine|fc1:"));
    }

    [Test]
    public void CanParseMatchingFragmentsAndMatches()
    {
        var exp = _dataSet.Experiments[0];
        var hit = exp.Analyses[0].Results.HitLists[0].Hits[0];
        var frag = hit.MatchingFragments.FirstOrDefault(f => f.Name == "B11");
        Assert.That(frag, Is.Not.Null);
        Assert.That(frag.TheoreticalMass, Is.EqualTo(1424.73869).Within(1e-6));
        Assert.That(frag.Matches, Has.Count.EqualTo(1));
        var match = frag.Matches[0];
        Assert.That(match.Id, Is.EqualTo("267319"));
        Assert.That(match.MassDifferenceDa, Is.EqualTo(0.00830999999993765).Within(1e-12));
        Assert.That(match.MassDifferencePpm, Is.EqualTo(5.83264851180369).Within(1e-10));
    }

    [Test]
    public void DataSet_Properties_Are_Correct()
    {
        Assert.That(_dataSet, Is.Not.Null);
        Assert.That(_dataSet.Version, Is.EqualTo("1.1"));
        Assert.That(_dataSet.Experiments, Is.Not.Null.And.Count.EqualTo(1));
    }

    [Test]
    public void MsMsExperiment_Properties_Are_Correct()
    {
        var exp = _dataSet.Experiments[0];
        Assert.That(exp.Id, Is.EqualTo("955"));
        Assert.That(exp.Source, Is.EqualTo(""));
        Assert.That(exp.Comment, Does.Contain("Characterization score: 221.13"));
    }

    [Test]
    public void InstrumentData_Properties_Are_Correct()
    {
        var exp = _dataSet.Experiments[0];
        var inst = exp.InstrumentData;
        Assert.That(inst, Is.Not.Null);
        Assert.That(inst.FragmentationMethod, Is.EqualTo("HCD"));
        Assert.That(inst.IonType, Is.EqualTo("BY"));
        Assert.That(inst.Intacts, Is.Not.Null.And.Count.EqualTo(1));
        Assert.That(inst.Fragments, Is.Not.Null.And.Not.Empty);

        var intact = inst.Intacts[0];
        Assert.That(intact.Id, Is.EqualTo("267562"));
        Assert.That(intact.MzMonoisotopic, Is.EqualTo(0));
        Assert.That(intact.MzAverage, Is.EqualTo(0));
        Assert.That(intact.MassMonoisotopic, Is.EqualTo(16153.783).Within(1e-6));
        Assert.That(intact.MassAverage, Is.EqualTo(0));
        Assert.That(intact.Intensity, Is.EqualTo(88912.86).Within(1e-2));

        var frag = inst.Fragments.First(f => f.Id == "267290");
        Assert.That(frag.MzMonoisotopic, Is.EqualTo(1795.87524).Within(1e-6));
        Assert.That(frag.MzAverage, Is.EqualTo(0));
        Assert.That(frag.MassMonoisotopic, Is.EqualTo(16153.8117).Within(1e-6));
        Assert.That(frag.MassAverage, Is.EqualTo(0));
        Assert.That(frag.Intensity, Is.EqualTo(86691.19).Within(1e-2));
    }

    [Test]
    public void Analysis_And_SearchParameters_Are_Correct()
    {
        var exp = _dataSet.Experiments[0];
        Assert.That(exp.Analyses, Is.Not.Null.And.Count.EqualTo(1));
        var analysis = exp.Analyses[0];
        Assert.That(analysis.Id, Is.EqualTo("1"));
        Assert.That(analysis.Type, Is.EqualTo("absolute_mass"));

        var search = analysis.SearchParameters;
        Assert.That(search, Is.Not.Null);
        Assert.That(search.Database, Is.Not.Null);
        Assert.That(search.Database.InternalName, Is.EqualTo("pseudomonas_aeruginosa_2012_11_top_down_complex"));
        Assert.That(search.Database.Display, Is.EqualTo("Pseudomonas Aeruginosa 2012 11 Top Down Complex"));
        Assert.That(search.PtmList, Does.Contain("Pyroglutamic acid from glutamic acid"));
        Assert.That(search.FixedModificationList, Is.Not.Null.And.Empty);
        Assert.That(search.TerminalModificationList, Is.Not.Null.And.Empty);
        Assert.That(search.IntactTolerance, Is.EqualTo("1000"));
        Assert.That(search.IntactToleranceUnit, Is.EqualTo("Da"));
        Assert.That(search.IntactMassType, Is.EqualTo("monoisotopic"));
        Assert.That(search.FragmentTolerance, Is.EqualTo("15"));
        Assert.That(search.FragmentToleranceUnit, Is.EqualTo("ppm"));
        Assert.That(search.FragmentMassType, Is.EqualTo("monoisotopic"));
        Assert.That(search.DeltaM, Is.True);
        Assert.That(search.Multiplexing, Is.False);
        Assert.That(search.ITraq, Is.False);
        Assert.That(search.Disulfide, Is.False);

        var hitCriteria = search.HitCriteria;
        Assert.That(hitCriteria, Is.Not.Null);
        Assert.That(hitCriteria.Max, Is.EqualTo(25));
        Assert.That(hitCriteria.MinimumMatchesNum, Is.EqualTo(4));
        Assert.That(hitCriteria.MinimumMatchesPercent, Is.EqualTo(0));
        Assert.That(hitCriteria.Score, Is.EqualTo("0"));
    }

    [Test]
    public void Results_And_HitList_Are_Correct()
    {
        var exp = _dataSet.Experiments[0];
        var results = exp.Analyses[0].Results;
        Assert.That(results, Is.Not.Null);
        Assert.That(results.HitLists, Is.Not.Null.And.Count.EqualTo(1));
        Assert.That(results.Legend, Is.EqualTo("Cysteine|fc1:"));

        var hitList = results.HitLists[0];
        Assert.That(hitList.IntactId, Is.EqualTo("1"));
        Assert.That(hitList.Hits, Is.Not.Null.And.Count.EqualTo(1));
    }

    [Test]
    public void Identification_And_MatchingSequence_Are_Correct()
    {
        var hit = _dataSet.Experiments[0].Analyses[0].Results.HitLists[0].Hits[0];
        Assert.That(hit.Id, Is.EqualTo("1"));
        Assert.That(hit.MatchingGeneId, Is.EqualTo("0"));
        Assert.That(hit.ProteinForm, Is.EqualTo("0"));
        Assert.That(hit.Description, Does.Contain("30S ribosomal protein S6"));
        Assert.That(hit.Signalp, Is.False);
        Assert.That(hit.Propep, Is.False);
        Assert.That(hit.SequenceLength, Is.EqualTo(139));
        Assert.That(hit.TheoreticalMass, Is.EqualTo(16154.70955469).Within(1e-8));
        Assert.That(hit.MassDifferenceDa, Is.EqualTo(-0.926554689971454).Within(1e-10));
        Assert.That(hit.MassDifferencePpm, Is.EqualTo(-57.3550819242343).Within(1e-10));
        Assert.That(hit.Score, Is.Not.Null);

        var seq = hit.MatchingSequence;
        Assert.That(seq, Is.Not.Null);
        Assert.That(seq.Resid, Does.StartWith("MRHYEIVFLVHPDQSEQV"));
        Assert.That(seq.Format, Does.Contain("M~RHYEIVFLVH|PDQSEQVGGMVERYTKAIEEDGGKIHRLE|D|WGRRQLAYAINNVHK|AH|YVLMNVECSAKALAELED|NFRYN|D|AVIRNLVMRRD|EAVTEQSEMLKAEESRNERRERRERPNDNAEGADGDDNSD|SDNADE"));
        Assert.That(seq.Display, Does.Contain("0@000M!0^000R!0@000H!0@000Y!0@000E!0@000I!0@000V!0@000F!0@000L!0@000V!0@000H%0@000P!0@000D!0@000Q!0@000S!0@000E!0@000Q!0@000V!0@000G!0@000G!0@000M!0@000V!0@000E!0@000R!0@000Y!0@000T!0@000K!0@000A!0@000I!0@000E!0@000E!0@000D!0@000G!0@000G!0@000K!0@000I!0@000H!0@000R!0@000L!0@000E%0@000D%0@000W!0@000G!0@000R!0@000R!0@000Q!0@000L!0@000A!0@000Y!0@000A!0@000I!0@000N!0@000N!0@000V!0@000H!0@000K%0@000A!0@000H%0@000Y!0@000V!0@000L!0@000M!0@000N!0@000V!0@000E!0@fc1C!0@000S!0@000A!0@000K!0@000A!0@000L!0@000A!0@000E!0@000L!0@000E!0@000D%0@000N!0@000F!0@000R!0@000Y!0@000N%0@000D%0@000A!0@000V!0@000I!0@000R!0@000N!0@000L!0@000V!0@000M!0@000R!0@000R!0@000D%0@000E!0@000A!0@000V!0@000T!0@000E!0@000Q!0@000S!0@000E!0@000M!0@000L!0@000K!0@000A!0@000E!0@000E!0@000S!0@000R!0@000N!0@000E!0@000R!0@000R!0@000E!0@000R!0@000R!0@000E!0@000R!0@000P!0@000N!0@000D!0@000N!0@000A!0@000E!0@000G!0@000A!0@000D!0@000G!0@000D!0@000D!0@000N!0@000S!0@000D%0@000S!0@000D!0@000N!0@000A!0@000D!0@000E!"));
    }

    [Test]
    public void Score_And_MatchingFragments_Are_Correct()
    {
        var hit = _dataSet.Experiments[0].Analyses[0].Results.HitLists[0].Hits[0];
        var score = hit.Score;
        Assert.That(score.PScore, Is.EqualTo("4.99e-012"));
        Assert.That(score.Expected, Is.EqualTo("1.51e-007"));
        Assert.That(score.McluckeyScore, Is.EqualTo("0.0"));

        Assert.That(hit.MatchingFragments, Is.Not.Null.And.Not.Empty);
        var frag = hit.MatchingFragments.First(f => f.Name == "B11");
        Assert.That(frag.Name, Is.EqualTo("B11"));
        Assert.That(frag.TheoreticalMass, Is.EqualTo(1424.73869).Within(1e-6));
        Assert.That(frag.Matches, Is.Not.Null.And.Count.EqualTo(1));
        var match = frag.Matches[0];
        Assert.That(match.Id, Is.EqualTo("267319"));
        Assert.That(match.MassDifferenceDa, Is.EqualTo(0.00830999999993765).Within(1e-12));
        Assert.That(match.MassDifferencePpm, Is.EqualTo(5.83264851180369).Within(1e-10));
    }
}
