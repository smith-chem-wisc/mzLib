using System.Text.RegularExpressions;

namespace UsefulProteomicsDatabases
{
    public class FastaHeaderFieldRegex
    {
        public FastaHeaderFieldRegex(string fieldName, string regularExpression, int match, int group)
        {
            FieldName = fieldName;
            Regex = new Regex(regularExpression);
            Match = match;
            Group = group;
        }

        public string FieldName { get; }

        public Regex Regex { get; }

        public int Match { get; }

        public int Group { get; }
    }
}