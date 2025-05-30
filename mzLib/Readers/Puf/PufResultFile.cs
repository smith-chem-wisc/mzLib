using System.Globalization;
using System.Xml.Linq;

namespace Readers.Puf;

/// <summary>
/// ResultFile implementation for Prosight PUF files.
/// Each result is a PufIdentification (i.e., a single identification from the file).
/// </summary>
internal class PufResultFile : ResultFile<PufMsMsExperiment>
{
    internal PufDataSet DataSet;

    public override SupportedFileType FileType => SupportedFileType.Puf;

    private Software _software = Software.ProsightPC;
    public override Software Software
    {
        get => _software;
        set => _software = value;
    }

    public PufResultFile() : base("") { }

    public PufResultFile(string filePath) : base(filePath, Software.ProsightPD) { }

    /// <summary>
    /// Loads all identifications from the PUF file into the Results list.
    /// </summary>
    public override void LoadResults()
    {
        if (!File.Exists(FilePath))
            throw new FileNotFoundException($"PUF file not found: {FilePath}");

        var stream = File.OpenRead(FilePath);
        DataSet = Read(stream);
        Results = DataSet.Experiments;
    }

    /// <summary>
    /// Writes the results to a specified output path as a simple tab-delimited file.
    /// </summary>
    public override void WriteResults(string outputPath)
    {
        var stream = File.Create(outputPath);
        Write(DataSet, stream);
    }

    internal static PufDataSet Read(Stream stream)
    {
        var doc = XDocument.Load(stream);
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

    internal static void Write(PufDataSet dataSet, Stream stream)
    {
        var doc = new XDocument(
            new XElement("data_set",
                new XAttribute("version", dataSet.Version ?? "1.1"),
                dataSet.Experiments.Select(exp =>
                    new XElement("ms-ms_experiment",
                        new XAttribute("id", exp.Id ?? ""),
                        new XAttribute("source", exp.Source ?? ""),
                        new XElement("comment", exp.Comment ?? ""),
                        exp.InstrumentData == null ? null :
                        new XElement("instrument_data",
                            new XElement("fragmentation_method", exp.InstrumentData.FragmentationMethod ?? ""),
                            new XElement("ion_type", exp.InstrumentData.IonType ?? ""),
                            new XElement("intact_list",
                                exp.InstrumentData.Intacts.Select(intact =>
                                    new XElement("intact",
                                        new XAttribute("id", intact.Id ?? ""),
                                        new XElement("mz_monoisotopic", intact.MzMonoisotopic ?? 0),
                                        new XElement("mz_average", intact.MzAverage ?? 0),
                                        new XElement("mass_monoisotopic", intact.MassMonoisotopic ?? 0),
                                        new XElement("mass_average", intact.MassAverage ?? 0),
                                        new XElement("intensity", intact.Intensity ?? 0)
                                    )
                                )
                            ),
                            new XElement("fragment_list",
                                exp.InstrumentData.Fragments.Select(frag =>
                                    new XElement("fragment",
                                        new XAttribute("id", frag.Id ?? ""),
                                        new XElement("mz_monoisotopic", frag.MzMonoisotopic ?? 0),
                                        new XElement("mz_average", frag.MzAverage ?? 0),
                                        new XElement("mass_monoisotopic", frag.MassMonoisotopic ?? 0),
                                        new XElement("mass_average", frag.MassAverage ?? 0),
                                        new XElement("intensity", frag.Intensity ?? 0)
                                    )
                                )
                            )
                        ),
                        exp.Analyses.Select(analysis =>
                            new XElement("analysis",
                                new XAttribute("id", analysis.Id ?? ""),
                                new XAttribute("type", analysis.Type ?? ""),
                                WriteSearchParameters(analysis.SearchParameters),
                                WriteResults(analysis.Results)
                            )
                        )
                    )
                )
            )
        );
        doc.Save(stream);
    }

    // --- Helper methods for parsing and writing sub-objects ---
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

    private static XElement WriteSearchParameters(PufSearchParameters search)
    {
        if (search == null) return null;
        return new XElement("search_parameters",
            search.Database == null ? null :
                new XElement("database",
                    new XElement("internal_name", search.Database.InternalName ?? ""),
                    new XElement("display", search.Database.Display ?? "")
                ),
            new XElement("ptm_list", search.PtmList.Select(ptm => new XElement("ptm", ptm))),
            new XElement("fixed_modification_list", search.FixedModificationList.Select(fm => new XElement("mod", fm))),
            new XElement("terminal_modification_list", search.TerminalModificationList.Select(tm => new XElement("mod", tm))),
            new XElement("intact_tolerance", new XAttribute("unit", search.IntactToleranceUnit ?? ""), search.IntactTolerance ?? ""),
            new XElement("intact_mass_type", search.IntactMassType ?? ""),
            new XElement("fragment_tolerance", new XAttribute("unit", search.FragmentToleranceUnit ?? ""), search.FragmentTolerance ?? ""),
            new XElement("fragment_mass_type", search.FragmentMassType ?? ""),
            new XElement("delta_m", search.DeltaM?.ToString().ToLowerInvariant() ?? ""),
            new XElement("multiplexing", search.Multiplexing?.ToString().ToLowerInvariant() ?? ""),
            new XElement("iTraq", search.ITraq?.ToString().ToLowerInvariant() ?? ""),
            new XElement("disulfide", search.Disulfide?.ToString().ToLowerInvariant() ?? ""),
            search.HitCriteria == null ? null :
                new XElement("hit_criteria",
                    new XAttribute("max", search.HitCriteria.Max?.ToString() ?? ""),
                    new XElement("minimum_matches_num", search.HitCriteria.MinimumMatchesNum?.ToString() ?? ""),
                    new XElement("minimum_matches_percent", search.HitCriteria.MinimumMatchesPercent?.ToString() ?? ""),
                    new XElement("score", search.HitCriteria.Score ?? "")
                )
        );
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

    private static XElement WriteResults(PufResults results)
    {
        if (results == null) return null;
        return new XElement("results",
            results.HitLists.Select(hitList =>
                new XElement("hit_list",
                    new XAttribute("intact_id", hitList.IntactId ?? ""),
                    hitList.Hits.Select(hit =>
                        new XElement("hit",
                            new XAttribute("id", hit.Id ?? ""),
                            new XElement("matching_gene_id", hit.MatchingGeneId ?? ""),
                            new XElement("protein_form", hit.ProteinForm ?? ""),
                            new XElement("description", hit.Description ?? ""),
                            new XElement("signalp", hit.Signalp?.ToString().ToLowerInvariant() ?? ""),
                            new XElement("propep", hit.Propep?.ToString().ToLowerInvariant() ?? ""),
                            hit.MatchingSequence == null ? null :
                                new XElement("matching_sequence",
                                    new XElement("resid", hit.MatchingSequence.Resid ?? ""),
                                    new XElement("format", hit.MatchingSequence.Format ?? ""),
                                    new XElement("display", hit.MatchingSequence.Display ?? "")
                                ),
                            new XElement("sequence_length", hit.SequenceLength?.ToString() ?? ""),
                            new XElement("theoretical_mass", hit.TheoreticalMass?.ToString(CultureInfo.InvariantCulture) ?? ""),
                            new XElement("mass_difference_da", hit.MassDifferenceDa?.ToString(CultureInfo.InvariantCulture) ?? ""),
                            new XElement("mass_difference_ppm", hit.MassDifferencePpm?.ToString(CultureInfo.InvariantCulture) ?? ""),
                            hit.Score == null ? null :
                                new XElement("score",
                                    new XElement("p_score", hit.Score.PScore ?? ""),
                                    new XElement("expected", hit.Score.Expected ?? ""),
                                    new XElement("mcluckey_score", hit.Score.McluckeyScore ?? "")
                                ),
                            new XElement("matching_fragment_list",
                                hit.MatchingFragments.Select(mf =>
                                    new XElement("matching_fragment",
                                        new XElement("name", mf.Name ?? ""),
                                        new XElement("theoretical_mass", mf.TheoreticalMass?.ToString(CultureInfo.InvariantCulture) ?? ""),
                                        mf.Matches.Select(match =>
                                            new XElement("match",
                                                new XAttribute("id", match.Id ?? ""),
                                                new XElement("mass_difference_da", match.MassDifferenceDa?.ToString(CultureInfo.InvariantCulture) ?? ""),
                                                new XElement("mass_difference_ppm", match.MassDifferencePpm?.ToString(CultureInfo.InvariantCulture) ?? "")
                                            )
                                        )
                                    )
                                )
                            )
                        )
                    )
                )
            ),
            new XElement("legend", results.Legend ?? "")
        );
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
