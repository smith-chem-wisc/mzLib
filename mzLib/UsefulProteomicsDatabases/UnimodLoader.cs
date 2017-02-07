using Chemistry;
using System;
using System.Collections.Generic;
using System.IO;
using System.Xml.Serialization;

namespace UsefulProteomicsDatabases
{
    internal class UnimodLoader
    {
        private static readonly Dictionary<string, string> DictOfElements = new Dictionary<string, string>
        {
            {"2H", "H{2}" },
            {"13C", "C{13}" },
            {"18O", "O{18}" },
            {"15N", "N{15}" }
        };

        internal static IEnumerable<Modification> ReadMods(string unimodLocation)
        {
            var unimodSerializer = new XmlSerializer(typeof(Generated.unimod_t));
            var deserialized = unimodSerializer.Deserialize(new FileStream(unimodLocation, FileMode.Open)) as Generated.unimod_t;

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
                        string val;
                        if (DictOfElements.TryGetValue(el.symbol, out val))
                        {
                            var tempCF = ChemicalFormula.ParseFormula(val);
                            tempCF.Multiply(int.Parse(el.number));
                            cf.Add(tempCF);
                        }
                        else
                        {
                            cf = null;
                            break;
                        }
                    }
                }

                var mm = cool.delta.mono_mass;

                foreach (var nice in cool.specificity)
                {
                    var tg = nice.site;
                    var pos = nice.position;
                    if (nice.NeutralLoss == null)
                        yield return new ModificationWithMassAndCf(id, new Tuple<string, string>("unimod", ac.ToString()), tg, pos, cf, mm, 0);
                    else
                        foreach (var nl in nice.NeutralLoss)
                            yield return new ModificationWithMassAndCf(id, new Tuple<string, string>("unimod", ac.ToString()), tg, pos, cf, mm, nl.mono_mass);

                }
            }
        }
    }
}