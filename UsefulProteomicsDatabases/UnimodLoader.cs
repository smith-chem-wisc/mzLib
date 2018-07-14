using Chemistry;
using Proteomics;
using System.Collections.Generic;
using System.IO;
using System.Xml.Serialization;
using UsefulProteomicsDatabases.Generated;

namespace UsefulProteomicsDatabases
{
    internal static class UnimodLoader
    {
        private static readonly Dictionary<string, string> DictOfElements = new Dictionary<string, string>
        {
            {"2H", "H{2}" },
            {"13C", "C{13}" },
            {"18O", "O{18}" },
            {"15N", "N{15}" }
        };

        private static readonly Dictionary<position_t, TerminusLocalization> positionDict = new Dictionary<position_t, TerminusLocalization>
            {
            {position_t.AnyCterm, TerminusLocalization.PepC },
            {position_t.ProteinCterm, TerminusLocalization.ProtC },
            {position_t.Anywhere, TerminusLocalization.Any },
            {position_t.AnyNterm, TerminusLocalization.NPep },
            {position_t.ProteinNterm, TerminusLocalization.NProt }
            };

        internal static IEnumerable<ModificationGeneral> ReadMods(string unimodLocation)
        {
            var unimodSerializer = new XmlSerializer(typeof(unimod_t));
            var deserialized = unimodSerializer.Deserialize(new FileStream(unimodLocation, FileMode.Open, FileAccess.Read, FileShare.Read)) as unimod_t;

            Dictionary<string, string> positionConversion = new Dictionary<string, string>()
            {
                { "Anywhere", "Anywhere."},
                { "AnyNterm", "Peptide N-terminal."},
                { "AnyCterm", "Peptide C-terminal."},
                { "ProteinNterm", "N-terminal."},
                { "ProteinCterm", "C-terminal."}
            };

            foreach (var cool in deserialized.modifications)
            {
                var id = cool.title;
                var ac = cool.record_id;
                ChemicalFormula cf = new ChemicalFormula();
                foreach (var el in cool.delta.element)
                {
                    try
                    {
                        cf.Add(el.symbol, int.Parse(el.number));
                    }
                    catch
                    {
                        var tempCF = ChemicalFormula.ParseFormula(DictOfElements[el.symbol]);
                        tempCF.Multiply(int.Parse(el.number));
                        cf.Add(tempCF);
                    }
                }

                foreach (var nice in cool.specificity)
                {
                    var tg = nice.site;
                    if (tg.Length > 1)
                        tg = "X"; //I think that we should allow motifs here using the trygetmotif
                    ModificationMotif.TryGetMotif(tg, out ModificationMotif motif);

                    if (positionConversion.TryGetValue(nice.position.ToString(), out string pos))
                    {
                        //do nothing, the new string value should be there
                    }
                    else
                    {
                        pos = null;
                    }

                    Dictionary<string, IList<string>> dblinks = new Dictionary<string, IList<string>>
                    {
                        { "Unimod",  new List<string>{ac.ToString() } },
                    };

                    if (nice.NeutralLoss == null)
                        yield return new ModificationGeneral(_Id: id, _ModificationType: "Unimod", _Target: motif, _Position: pos, _ChemicalFormula: cf, _DatabaseReference: dblinks);
                    else
                    {
                        Dictionary<MassSpectrometry.DissociationType, List<double>> neutralLosses = null;
                        foreach (var nl in nice.NeutralLoss)
                        {
                            ChemicalFormula cfnl = new ChemicalFormula();
                            if (nl.mono_mass == 0)
                            {
                                if (neutralLosses == null)
                                {
                                    neutralLosses = new Dictionary<MassSpectrometry.DissociationType, List<double>>();
                                    neutralLosses.Add(MassSpectrometry.DissociationType.AnyActivationType, new List<double> { 0 });
                                }
                                else
                                {
                                    if (neutralLosses.ContainsKey(MassSpectrometry.DissociationType.AnyActivationType))
                                    {
                                        if (!neutralLosses[MassSpectrometry.DissociationType.AnyActivationType].Contains(0))
                                        {
                                            neutralLosses[MassSpectrometry.DissociationType.AnyActivationType].Add(0);
                                        }
                                    }//we don't need an else cuz it's already there
                                }
                            }
                            else
                            {
                                foreach (var el in nl.element)
                                {
                                    try
                                    {
                                        cfnl.Add(el.symbol, int.Parse(el.number));
                                    }
                                    catch
                                    {
                                        var tempCF = ChemicalFormula.ParseFormula(DictOfElements[el.symbol]);
                                        tempCF.Multiply(int.Parse(el.number));
                                        cfnl.Add(tempCF);
                                    }
                                }
                                if (neutralLosses == null)
                                {
                                    neutralLosses = new Dictionary<MassSpectrometry.DissociationType, List<double>>();
                                    neutralLosses.Add(MassSpectrometry.DissociationType.AnyActivationType, new List<double> { cfnl.MonoisotopicMass });
                                }
                                else
                                {
                                    if (neutralLosses.ContainsKey(MassSpectrometry.DissociationType.AnyActivationType))
                                    {
                                        if (!neutralLosses[MassSpectrometry.DissociationType.AnyActivationType].Contains(cfnl.MonoisotopicMass))
                                        {
                                            neutralLosses[MassSpectrometry.DissociationType.AnyActivationType].Add(cfnl.MonoisotopicMass);
                                        }
                                    }//we don't need an else cuz it's already there
                                }
                            }
                        }
                        yield return new ModificationGeneral(_Id: id, _Target: motif, _Position: "Anywhere.", _ModificationType: "Unimod", _ChemicalFormula: cf, _DatabaseReference: dblinks, _NeutralLosses: neutralLosses);
                    }
                }
            }
        }
    }
}