using MzLibUtil;

namespace Readers
{
    /// <summary>
    /// The first line of MetaMorpheus's DIA tables: whether a contaminant database was searched, so that "no
    /// contaminants" is never true by default. A table without it is refused, not read as "not assessed".
    /// </summary>
    internal static class DiaContaminantsLine
    {
        private const string Prefix = "# contaminants: ";

        public static void Write(TextWriter writer, bool contaminantsAssessed) =>
            writer.WriteLine(Prefix + (contaminantsAssessed ? "assessed" : "not assessed"));

        /// <exception cref="MzLibException">The first line does not say whether contaminants were assessed.</exception>
        public static bool Read(TextReader reader, string filePath, string table) => (reader.ReadLine() ?? "") switch
        {
            Prefix + "assessed" => true,
            Prefix + "not assessed" => false,
            _ => throw new MzLibException($"{filePath} is not a {table}: its first line must say whether contaminants were assessed."),
        };
    }
}
