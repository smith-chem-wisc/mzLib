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

        internal static IEnumerable<ModificationWithLocation> ReadMods(string unimodLocation)
        {
            var unimodSerializer = new XmlSerializer(typeof(unimod_t));
            var deserialized = unimodSerializer.Deserialize(new FileStream(unimodLocation, FileMode.Open, FileAccess.Read, FileShare.Read)) as unimod_t;

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
                        tg = "X";
                    ModificationMotif.TryGetMotif(tg, out ModificationMotif motif);
                    var pos = nice.position;
                    IDictionary<string, IList<string>> dblinks = new Dictionary<string, IList<string>>
                    {
                        { "Unimod",  new List<string>{ac.ToString() } },
                    };

                    if (nice.NeutralLoss == null)
                        yield return new ModificationWithMassAndCf(id + " on " + motif + " at " + positionDict[pos], "Unimod", motif, positionDict[pos], cf, linksToOtherDbs: dblinks);
                    else
                    {
                        List<double> neutralLosses = new List<double>();
                        foreach (var nl in nice.NeutralLoss)
                        {
                            ChemicalFormula cfnl = new ChemicalFormula();
                            if (nl.mono_mass == 0)
                                neutralLosses.Add(0);
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
                                neutralLosses.Add(cfnl.MonoisotopicMass);
                            }
                        }
                        yield return new ModificationWithMassAndCf(id + " on " + motif + " at " + positionDict[pos], "Unimod", motif, positionDict[pos], cf, linksToOtherDbs: dblinks, neutralLosses: neutralLosses);
                    }
                }
            }
        }
    }
}