using System.Globalization;
using System.Xml.Linq;

namespace Readers.Puf;

public static class PufParser
{
    internal static PufDataSet Parse(string filePath)
    {
        var doc = XDocument.Load(filePath);
        var dataSetElem = doc.Element("data_set");
        var dataSet = new PufDataSet
        {
            Version = dataSetElem?.Attribute("version")?.Value
        };

        foreach (var expElem in dataSetElem.Elements("ms-ms_experiment"))
        {
            var exp = new PufMsMsExperiment
            {
                Id = expElem.Attribute("id")?.Value,
                Source = expElem.Attribute("source")?.Value,
                Comment = expElem.Element("comment")?.Value
            };

            var instElem = expElem.Element("instrument_data");
            if (instElem != null)
            {
                var inst = new PufDataScan
                {
                    FragmentationMethod = instElem.Element("fragmentation_method")?.Value,
                    IonType = instElem.Element("ion_type")?.Value
                };

                var intactListElem = instElem.Element("intact_list");
                if (intactListElem != null)
                {
                    foreach (var intactElem in intactListElem.Elements("intact"))
                    {
                        inst.Intacts.Add(new PufIntact
                        {
                            Id = intactElem.Attribute("id")?.Value,
                            MzMonoisotopic = ParseNullableDouble(intactElem.Element("mz_monoisotopic")?.Value),
                            MzAverage = ParseNullableDouble(intactElem.Element("mz_average")?.Value),
                            MassMonoisotopic = ParseNullableDouble(intactElem.Element("mass_monoisotopic")?.Value),
                            MassAverage = ParseNullableDouble(intactElem.Element("mass_average")?.Value),
                            Intensity = ParseNullableDouble(intactElem.Element("intensity")?.Value)
                        });
                    }
                }

                var fragListElem = instElem.Element("fragment_list");
                if (fragListElem != null)
                {
                    foreach (var fragElem in fragListElem.Elements("fragment"))
                    {
                        inst.Fragments.Add(new PufFragment
                        {
                            Id = fragElem.Attribute("id")?.Value,
                            MzMonoisotopic = ParseNullableDouble(fragElem.Element("mz_monoisotopic")?.Value),
                            MzAverage = ParseNullableDouble(fragElem.Element("mz_average")?.Value),
                            MassMonoisotopic = ParseNullableDouble(fragElem.Element("mass_monoisotopic")?.Value),
                            MassAverage = ParseNullableDouble(fragElem.Element("mass_average")?.Value),
                            Intensity = ParseNullableDouble(fragElem.Element("intensity")?.Value)
                        });
                    }
                }

                exp.InstrumentData = inst;
            }

            foreach (var analysisElem in expElem.Elements("analysis"))
            {
                var analysis = new PufAnalysis
                {
                    Id = analysisElem.Attribute("id")?.Value,
                    Type = analysisElem.Attribute("type")?.Value,
                    SearchParameters = ParseSearchParameters(analysisElem.Element("search_parameters")),
                    Results = ParseResults(analysisElem.Element("results"))
                };

                exp.Analyses.Add(analysis);
            }

            dataSet.Experiments.Add(exp);
        }

        return dataSet;
    }

    private static PufSearchParameters ParseSearchParameters(XElement elem)
    {
        if (elem == null) return null;
        var dbElem = elem.Element("database");
        var ptmListElem = elem.Element("ptm_list");
        var fixedModElem = elem.Element("fixed_modification_list");
        var terminalModElem = elem.Element("terminal_modification_list");
        var hitCriteriaElem = elem.Element("hit_criteria");

        return new PufSearchParameters
        {
            Database = dbElem == null ? null : new PufDatabase
            {
                InternalName = dbElem.Element("internal_name")?.Value,
                Display = dbElem.Element("display")?.Value
            },
            PtmList = ptmListElem?.Elements().Select(e => e.Value).ToList() ?? new List<string>(),
            FixedModificationList = fixedModElem?.Elements().Select(e => e.Value).ToList() ?? new List<string>(),
            TerminalModificationList = terminalModElem?.Elements().Select(e => e.Value).ToList() ?? new List<string>(),
            IntactTolerance = elem.Element("intact_tolerance")?.Value,
            IntactToleranceUnit = elem.Element("intact_tolerance")?.Attribute("unit")?.Value,
            IntactMassType = elem.Element("intact_mass_type")?.Value,
            FragmentTolerance = elem.Element("fragment_tolerance")?.Value,
            FragmentToleranceUnit = elem.Element("fragment_tolerance")?.Attribute("unit")?.Value,
            FragmentMassType = elem.Element("fragment_mass_type")?.Value,
            DeltaM = ParseNullableBool(elem.Element("delta_m")?.Value),
            Multiplexing = ParseNullableBool(elem.Element("multiplexing")?.Value),
            ITraq = ParseNullableBool(elem.Element("iTraq")?.Value),
            Disulfide = ParseNullableBool(elem.Element("disulfide")?.Value),
            HitCriteria = hitCriteriaElem == null ? null : new PufHitCriteria
            {
                Max = ParseNullableInt(hitCriteriaElem.Attribute("max")?.Value),
                MinimumMatchesNum = ParseNullableInt(hitCriteriaElem.Element("minimum_matches_num")?.Value),
                MinimumMatchesPercent = ParseNullableInt(hitCriteriaElem.Element("minimum_matches_percent")?.Value),
                Score = hitCriteriaElem.Element("score")?.Value
            }
        };
    }

    private static PufResults ParseResults(XElement elem)
    {
        if (elem == null) return null;
        var results = new PufResults();
        foreach (var hitListElem in elem.Elements("hit_list"))
        {
            var hitList = new PufIdentificationCollection
            {
                IntactId = hitListElem.Attribute("intact_id")?.Value
            };
            foreach (var hitElem in hitListElem.Elements("hit"))
            {
                var hit = new PufIdentification
                {
                    Id = hitElem.Attribute("id")?.Value,
                    MatchingGeneId = hitElem.Element("matching_gene_id")?.Value,
                    ProteinForm = hitElem.Element("protein_form")?.Value,
                    Description = hitElem.Element("description")?.Value,
                    Signalp = ParseNullableBool(hitElem.Element("signalp")?.Value),
                    Propep = ParseNullableBool(hitElem.Element("propep")?.Value),
                    MatchingSequence = hitElem.Element("matching_sequence") == null ? null : new PufMatchingSequence
                    {
                        Resid = hitElem.Element("matching_sequence").Element("resid")?.Value,
                        Format = hitElem.Element("matching_sequence").Element("format")?.Value,
                        Display = hitElem.Element("matching_sequence").Element("display")?.Value
                    },
                    SequenceLength = ParseNullableInt(hitElem.Element("sequence_length")?.Value),
                    TheoreticalMass = ParseNullableDouble(hitElem.Element("theoretical_mass")?.Value),
                    MassDifferenceDa = ParseNullableDouble(hitElem.Element("mass_difference_da")?.Value),
                    MassDifferencePpm = ParseNullableDouble(hitElem.Element("mass_difference_ppm")?.Value),
                    Score = hitElem.Element("score") == null ? null : new PufScore
                    {
                        PScore = hitElem.Element("score").Element("p_score")?.Value,
                        Expected = hitElem.Element("score").Element("expected")?.Value,
                        McluckeyScore = hitElem.Element("score").Element("mcluckey_score")?.Value
                    }
                };

                // Matching fragments
                var matchingFragmentListElem = hitElem.Element("matching_fragment_list");
                if (matchingFragmentListElem != null)
                {
                    foreach (var fragElem in matchingFragmentListElem.Elements("matching_fragment"))
                    {
                        var matchingFragment = new PufMatchingFragment
                        {
                            Name = fragElem.Element("name")?.Value,
                            TheoreticalMass = ParseNullableDouble(fragElem.Element("theoretical_mass")?.Value)
                        };
                        foreach (var matchElem in fragElem.Elements("match"))
                        {
                            matchingFragment.Matches.Add(new PufFragmentMatch
                            {
                                Id = matchElem.Attribute("id")?.Value,
                                MassDifferenceDa = ParseNullableDouble(matchElem.Element("mass_difference_da")?.Value),
                                MassDifferencePpm = ParseNullableDouble(matchElem.Element("mass_difference_ppm")?.Value)
                            });
                        }
                        hit.MatchingFragments.Add(matchingFragment);
                    }
                }

                hitList.Hits.Add(hit);
            }
            results.HitLists.Add(hitList);
        }
        results.Legend = elem.Element("legend")?.Value;
        return results;
    }

    private static bool? ParseNullableBool(string s)
    {
        if (string.IsNullOrWhiteSpace(s)) return null;
        if (bool.TryParse(s, out var b)) return b;
        if (s == "1" || s.Equals("true", StringComparison.OrdinalIgnoreCase) || s.Equals("yes", StringComparison.OrdinalIgnoreCase)) return true;
        if (s == "0" || s.Equals("false", StringComparison.OrdinalIgnoreCase) || s.Equals("no", StringComparison.OrdinalIgnoreCase)) return false;
        return null;
    }

    private static int? ParseNullableInt(string s)
    {
        if (int.TryParse(s, out var i)) return i;
        return null;
    }

    private static double? ParseNullableDouble(string s)
    {
        if (double.TryParse(s, NumberStyles.Any, CultureInfo.InvariantCulture, out var d)) return d;
        return null;
    }
}