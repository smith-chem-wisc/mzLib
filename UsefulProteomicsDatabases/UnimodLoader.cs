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

        private static readonly Dictionary<position_t, ModLocationOnPeptideOrProtein> positionDict = new Dictionary<position_t, ModLocationOnPeptideOrProtein>
            {
            {position_t.AnyCterm, ModLocationOnPeptideOrProtein.PepC },
            {position_t.ProteinCterm, ModLocationOnPeptideOrProtein.ProtC },
            {position_t.Anywhere, ModLocationOnPeptideOrProtein.Any },
            {position_t.AnyNterm, ModLocationOnPeptideOrProtein.NPep },
            {position_t.ProteinNterm, ModLocationOnPeptideOrProtein.NProt }
            };

        internal static IEnumerable<Modification> ReadMods(string unimodLocation)
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

            foreach (var mod in deserialized.modifications)
            {
                var id = mod.title;
                var ac = mod.record_id;
                ChemicalFormula cf = new ChemicalFormula();
                foreach (var el in mod.delta.element)
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

                //TODO Add "on motif" to the ID field
                foreach (var target in mod.specificity)
                {
                    var tg = target.site;
                    if (tg.Length > 1)
                        tg = "X"; //I think that we should allow motifs here using the trygetmotif
                    ModificationMotif.TryGetMotif(tg, out ModificationMotif motif);

                    if (positionConversion.TryGetValue(target.position.ToString(), out string pos))
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

                    if (target.NeutralLoss == null)
                        yield return new Modification(_originalId: id, _modificationType: "Unimod", _target: motif, _locationRestriction: pos, _chemicalFormula: cf, _databaseReference: dblinks);
                    else
                    {
                        Dictionary<MassSpectrometry.DissociationType, List<double>> neutralLosses = null;
                        foreach (var nl in target.NeutralLoss)
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
                        yield return new Modification(_originalId: id, _target: motif, _locationRestriction: "Anywhere.", _modificationType: "Unimod", _chemicalFormula: cf, _databaseReference: dblinks, _neutralLosses: neutralLosses);
                    }
                }
            }
        }
    }
}