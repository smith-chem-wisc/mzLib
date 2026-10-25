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

        /// <param name="displayedIsoformId">The isoform this sequence is (UniProt's &lt;isoform&gt; whose
        /// &lt;sequence type="displayed"/&gt;), or null (or empty) when the entry lists no isoforms.</param>
        public UniProtSequenceAttributes(int length, int mass, string checkSum, DateTime entryModified, int sequenceVersion, bool? isPrecursor, FragmentType fragment, string displayedIsoformId)
            : this(length, mass, checkSum, entryModified, sequenceVersion, isPrecursor, fragment)
        {
            DisplayedIsoformId = string.IsNullOrEmpty(displayedIsoformId) ? null : displayedIsoformId;
        }

        /// <summary>
        /// The UniProt isoform whose sequence this is, e.g. "P04406-1", read from the entry's
        /// &lt;comment type="alternative products"&gt;; null when the entry lists no isoforms (or came from FASTA).
        /// The displayed isoform can carry any number, so "-1" is never assumed. It decides whether a reference scoped
        /// to an isoform (<see cref="DatabaseReference.Molecule"/>) applies to this sequence.
        /// </summary>
        public string DisplayedIsoformId { get; private set; } //optional
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
