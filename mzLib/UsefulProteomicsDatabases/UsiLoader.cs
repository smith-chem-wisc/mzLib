using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace UsefulProteomicsDatabases
{
    public static class UsiLoader
    {
        public static bool TryGetSpectrum(string usi, out Object spectrum, UsiDatabase database = UsiDatabase.PRIDE)
        {
            string[] usiSplit = usi.Split(':');
            if (usiSplit.Length < 5 | usiSplit.Length > 6)
                throw new ArgumentException("USI string was improperly formatted");

            string baseUrl = GetRequestUrl(database);
            string apiQuery = baseUrl + "usi=" + string.Join('&', usiSplit);

            var apiResponse = Loaders.AwaitAsync_GetSomeData(apiQuery);

            spectrum = null;
            return false;
        }

        

        internal static string GetRequestUrl(UsiDatabase db)
        {
            switch (db)
            {
                case (UsiDatabase.PRIDE):
                    return "https://www.ebi.ac.uk/pride/ws/archive/v2/spectra?";
                default:
                    throw new NotImplementedException();
            }
        }
    }

    public enum UsiDatabase
    {
        PRIDE,
        Unknown
    }
}
