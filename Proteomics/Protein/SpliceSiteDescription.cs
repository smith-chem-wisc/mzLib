using System;
using System.Collections.Generic;
using System.Text;

namespace Proteomics
{
    public class SpliceSiteDescription
    {
        public SpliceSiteDescription(string description)
        {
            if (description == null)
            {
                return;
            }
            Description = description;

            // Parse description info
            string[] spliceFields = description.Split(@"\t");
            if (spliceFields.Length < 6) { return; }
            ChromosomeName = spliceFields[0];
            Strand = spliceFields[1];
            Start = int.Parse(spliceFields[2]);
            End = int.Parse(spliceFields[3]);
            NextStart = int.Parse(spliceFields[4]);
            NextEnd = int.Parse(spliceFields[5]);

            // Parse novel bool
            if (spliceFields.Length > 6)
            {
                _novel = bool.Parse(spliceFields[6]);
            }
        }

        public string Description { get; private set; }
        public string ChromosomeName { get; }
        public string Strand { get; }
        public int Start { get; }
        public int End { get; }
        public int NextStart { get; }
        public int NextEnd { get; }
        private bool? _novel;
        // Updates Description when setting Novel
        public bool? Novel { get { return _novel; } 
            set { 
                if (_novel == value) { return; }
                _novel = value;
                if (Description.Split(@"\t").Length == 7)
                {
                    Description = value == null ? Description.Substring(0, Description.LastIndexOf(@"\t", System.StringComparison.Ordinal)) :
                        Description.Substring(0, Description.LastIndexOf(@"\t", System.StringComparison.Ordinal)) + @"\t" + value;
                }
                else
                {
                    Description = Description + @"\t" + value;
                }
            } }

        /// <summary>
        /// Returns the description, may have been modified after construction by changing Novel 
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            return Description;
        }

        public override bool Equals(object obj)
        {
            SpliceSiteDescription s = obj as SpliceSiteDescription;
            return s != null && s.Description.Equals(Description);
        }

        public override int GetHashCode()
        {
            return (Description ?? "").GetHashCode();
        }
    }
}
