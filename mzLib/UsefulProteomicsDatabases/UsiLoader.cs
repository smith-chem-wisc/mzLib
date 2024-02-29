using Easy.Common.Extensions;
using MassSpectrometry;
using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Net.Http.Json;
using System.Xml.Linq;


namespace UsefulProteomicsDatabases
{
    public enum UsiDatabase
    {
        PRIDE, 
        PeptideAtlas,
        MassIVE,
        ProteomeExchange,
        JPOST
    }

    public static class UsiLoader
    {
        /// <summary>
        /// Returns a properly formatted api query for a given database
        /// As of 1/8/2024 : Proteome Exchange, PeptideAtlas, and MassIVE are working. PRIDE and JPOST only return 404, regardless of query.
        /// </summary>
        /// <exception cref="NotImplementedException"></exception>
        public static string GetApiQuery(string usi, UsiDatabase database = UsiDatabase.ProteomeExchange, string version = "0.1")
        {
            switch(database)
            {
                case (UsiDatabase.ProteomeExchange):
                    return $"http://proteomecentral.proteomexchange.org/api/proxi/v{version}/spectra?resultType=full&usi={usi}";
                case (UsiDatabase.PRIDE):
                    return $"https://www.ebi.ac.uk/pride/ws/archive/v{version}/spectra?resultType=full&usi={usi}";
                case (UsiDatabase.PeptideAtlas):
                    return $"http://www.peptideatlas.org/api/proxi/v{version}/spectra?resultType=full&usi={usi}";
                case (UsiDatabase.MassIVE):
                    return $"http://massive.ucsd.edu/ProteoSAFe/proxi/v{version}/spectra?resultType=full&usi={usi}";
                case (UsiDatabase.JPOST):
                    return $"https://repository.jpostdb.org/proxi/spectra?resultType=full&usi={usi}";
                default:
                    throw new NotImplementedException("Could not generate API query for the given database.");
            }
        }

        /// <summary>
        /// Fetches a list of spectra via USI
        /// </summary>
        /// <param name="usi"> Unique Spectrum Identifier </param>
        /// <param name="usiSpectraList"> Returns a list of all matching spectra as UsiSpectrum </param>
        /// <param name="database"> A website/service that hosts spectra. Default: ProteomeExchange </param>
        /// <returns> True if spectra was found and downloaded succesfully, false otherwise </returns>
        public static bool TryGetSpectra(string usi, out List<UsiSpectrumFromJSON> usiSpectraList, UsiDatabase database = UsiDatabase.ProteomeExchange)
        {
            usiSpectraList = null;
            string apiQuery = GetApiQuery(usi, database);
            var apiResponse = Loaders.AwaitAsync_GetSomeData(apiQuery).Result;
            if (!apiResponse.IsSuccessStatusCode)
                return false;
            string result = apiResponse.Content.ReadAsStringAsync().Result;
            usiSpectraList = JsonConvert.DeserializeObject<List<UsiSpectrumFromJSON>>(result);
            if (usiSpectraList.IsNotNullOrEmpty() && usiSpectraList.Any(spectra => spectra != null))
                return true;
            else
                return false;
        }

        /// <summary>
        /// Fetches a list of spectra via USI
        /// </summary>
        /// <param name="identifier"> The unique dataset identifier. ex: "PXD000561"</param>
        /// <param name="run"> Name of the file within the dataset </param>
        /// <param name="typeFlag"> Type of identifier. Most commonly "scan" </param>
        /// <param name="index"> Unique spectrum identifier. Most commonly scan number</param>
        /// <param name="interpretation">Optional spectrum iterpretation. ex: "VLHPLEGAVVIIFK/2"</param>
        /// <param name="usiSpectraList"> Returns a list of all matching spectra as UsiSpectrum </param>
        /// <param name="prefix"> Start of the usi. Default is "mzSpec" </param>
        /// <param name="database"> A website/service that hosts spectra. Default: ProteomeExchange </param>
        /// <returns> True if spectra was found and downloaded succesfully, false otherwise </returns>
        public static bool TryGetSpectra(string identifier, string run, string typeFlag, string index, string interpretation, 
            out List<UsiSpectrumFromJSON> usiSpectraList, string prefix = "mzspec", UsiDatabase database = UsiDatabase.ProteomeExchange)
        {
            string usi = string.Join(':', new List<string> { prefix, identifier, run, typeFlag, index, interpretation });
            return TryGetSpectra(usi, out usiSpectraList, database);
        }

        /// <summary>
        /// Fetches a list of spectra via USI
        /// </summary>
        /// <param name="identifier"> The unique dataset identifier. ex: "PXD000561"</param>
        /// <param name="run"> Name of the file within the dataset </param>
        /// <param name="typeFlag"> Type of identifier. Most commonly "scan" </param>
        /// <param name="index"> Unique spectrum identifier. Most commonly scan number</param>
        /// <param name="usiSpectraList"> Returns a list of all matching spectra as UsiSpectrum </param>
        /// <param name="prefix"> Start of the usi. Default is "mzSpec" </param>
        /// <param name="database"> A website/service that hosts spectra. Default: ProteomeExchange </param>
        /// <returns> True if spectra was found and downloaded succesfully, false otherwise </returns>
        public static bool TryGetSpectra(string identifier, string run, string typeFlag, string index, 
            out List<UsiSpectrumFromJSON> usiSpectraList, string prefix = "mzspec", UsiDatabase database = UsiDatabase.ProteomeExchange)
        {
            string usi = string.Join(':', new List<string> { prefix, identifier, run, typeFlag, index });
            return TryGetSpectra(usi, out usiSpectraList, database);
        }
    }

    public class UsiSpectrumFromJSON
    {
        [JsonProperty("attributes")]
        internal IList<UsiAttribute> Attributes { get; set; }
        [JsonProperty("intensities")]
        internal IList<double> Intensities { get; set; }
        [JsonProperty("mzs")]
        internal IList<double> Mzs { get; set; }

        public MzSpectrum GetMzSpectrum()
        {
            if (Intensities == null || Mzs == null)
                return null;
            if (Intensities.Count != Mzs.Count)
                throw new ArgumentException("Intensity array and mz array have unequal length");
            return new MzSpectrum(Mzs.ToArray(), Intensities.ToArray(), shouldCopy: true);
        }

        public double? GetIsolationWindowTargetMz()
        {
            var isolationAttribute = Attributes.Where(a => a.Accession.Equals("MS:1000827"));
            if (isolationAttribute.Count() != 1)
                return null;
            if (!double.TryParse(isolationAttribute.First().Value, out double windowTargetMz))
                return null;
            else
                return windowTargetMz;
        }

        public bool TryGetBaseSequence(out string sequence)
        {
            sequence = null;
            var sequenceAttribute = Attributes.Where(a => a.Accession.Equals("MS:1000888"));
            if (sequenceAttribute.Count() != 1)
                return false;
            sequence = sequenceAttribute.First().Value;
            return true;
        }
    }

    public class UsiAttribute
    {
        [JsonProperty("accession")]
        public string Accession { get; set; }
        [JsonProperty("name")]
        public string Name { get; set; }
        [JsonProperty("value")]
        public string Value { get; set; }

    }

}
