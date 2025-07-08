using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Proteomics
{
    public class UniProtSequenceAttributes
    {
        public int Length { get; set; } //mandatory
        public int Mass { get; set; } //mandatory
        public string CheckSum { get; set; } //mandatory
        public DateTime EntryModified { get; set; } //mandatory
        public int SequenceVersion { get; set; } //mandatory
        public bool? IsPrecursor { get; set; } //optional
        public FragmentType Fragment { get; set; } //optional
        public enum FragmentType
        {
            unspecified,
            single,
            multiple
        }
        public UniProtSequenceAttributes(int length, int mass, string checkSum, DateTime entryModified, int sequenceVersion, bool? isPrecursor = null, FragmentType fragment = FragmentType.unspecified)
        {
            Length = length;
            Mass = mass;
            CheckSum = checkSum;
            EntryModified = entryModified;
            SequenceVersion = sequenceVersion;
            IsPrecursor = isPrecursor;
            Fragment = fragment;
        }

    }
}
