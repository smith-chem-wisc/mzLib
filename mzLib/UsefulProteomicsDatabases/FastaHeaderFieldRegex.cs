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

        public string ApplyRegex(string input)
        {
            string? result = null;
            var matches = Regex.Matches(input);
            if (matches.Count > Match && matches[Match].Groups.Count > Group)
            {
                result = matches[Match].Groups[Group].Value;
            }

            return result!;
        }
    }
}