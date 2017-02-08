using Chemistry;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Xml.Serialization;
using UsefulProteomicsDatabases.Generated;

namespace UsefulProteomicsDatabases
{
    internal static class UnimodLoader
    {

        #region Private Fields

        private static readonly Dictionary<string, string> DictOfElements = new Dictionary<string, string>
        {
            {"2H", "H{2}" },
            {"13C", "C{13}" },
            {"18O", "O{18}" },
            {"15N", "N{15}" }
        };

        private static readonly Dictionary<position_t, ModificationSites> positionDict = new Dictionary<position_t, ModificationSites>
            {
            {position_t.AnyCterm, ModificationSites.PepC },
            {position_t.ProteinCterm, ModificationSites.ProtC },
            {position_t.Anywhere, ModificationSites.Any },
            {position_t.AnyNterm, ModificationSites.NPep },
            {position_t.ProteinNterm, ModificationSites.NProt }
            };

        #endregion Private Fields

        #region Internal Methods

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
                        var tempCF = ChemicalFormula.ParseFormula(DictOfElements[el.symbol]);
                        tempCF.Multiply(int.Parse(el.number));
                        cf.Add(tempCF);
                    }
                }

                var mm = cool.delta.mono_mass;

                foreach (var nice in cool.specificity)
                {
                    var tg = nice.site;
                    var pos = nice.position;
                    if (nice.NeutralLoss == null)
                        yield return new ModificationWithMassAndCf(id, new Tuple<string, string>("unimod", ac.ToString()), tg, positionDict[pos], cf, mm, 0, Path.GetFileNameWithoutExtension(unimodLocation));
                    else
                        foreach (var nl in nice.NeutralLoss)
                            yield return new ModificationWithMassAndCf(id, new Tuple<string, string>("unimod", ac.ToString()), tg, positionDict[pos], cf, mm, nl.mono_mass, Path.GetFileNameWithoutExtension(unimodLocation));
                }
            }
        }

        #endregion Internal Methods

    }
}