using Chemistry;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Proteomics
{
    public class UniProtSequenceAttributes
    {
        public int Length { get; private set; } //mandatory
        public int Mass { get; private set; } //mandatory
        public string Checksum { get; private set; } //mandatory
        public DateTime EntryModified { get; private set; } //mandatory
        public int SequenceVersion { get; private set; } //mandatory
        public bool? IsPrecursor { get; private set; } //optional
        public FragmentType Fragment { get; private set; } //optional
        public enum FragmentType
        {
            //note that these enum values do not follow the .Net naming convention (PascalCase) because they are specified by UniProt as all lowercase
            unspecified = 0,
            single = 1,
            multiple = 2
        }
        public UniProtSequenceAttributes(int length, int mass, string checkSum, DateTime entryModified, int sequenceVersion, bool? isPrecursor = null, FragmentType fragment = FragmentType.unspecified)
        {
            Length = length;
            Mass = mass;
            Checksum = checkSum;
            EntryModified = entryModified;
            SequenceVersion = sequenceVersion;
            IsPrecursor = isPrecursor;
            Fragment = fragment;
        }
        public void UpdateLengthAttribute(int newLength)
        {
            Length = newLength;
        }
        public void UpdateLengthAttribute(string newBaseSequence)
        {
            Length = newBaseSequence.Length;
        }
        public void UpdateMassAttribute(int newMass)
        {
            Mass = newMass;
        }
        public void UpdateMassAttribute(string newBaseSequence)
        {
            Mass = (int)Math.Round(new PeptideWithSetModifications(newBaseSequence, new Dictionary<string, Modification>()).MonoisotopicMass); 
        }
    }
}
