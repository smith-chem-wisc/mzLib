using System;
using Proteomics;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace UsefulProteomicsDatabases
{
    public class ProteinPeffEntry
    {
        public string Accession { get; set; }
        //public string Description { get; set; }
        //public string GeneName { get; set; }
        public string Organism { get; set; }
        public string Sequence { get; set; }
        //public List<string> Synonyms { get; set; } = new List<string>();
        // List<string> CrossReferences { get; set; } = new List<string>();
        //public List<string> Features { get; set; } = new List<string>();
        public PeffEntrySourceDatabase SourceDatabase { get; set; } = PeffEntrySourceDatabase.Uniprot;

        public ProteinPeffEntry(string accession, string description, string geneName, string organism, string sequence, PeffEntrySourceDatabase source=PeffEntrySourceDatabase.Uniprot)
        {
            Accession = accession;
            //Description = description;
            //GeneName = geneName;
            Organism = organism;
            Sequence = sequence;
            SourceDatabase = source;
        }
        public ProteinPeffEntry(ProteinXmlEntry proteinXmlEntry)
        {
            Accession = proteinXmlEntry.Accession;
            //Description = proteinXmlEntry.Description;
            //GeneName = proteinXmlEntry.GeneName;
            Organism = proteinXmlEntry.Organism;
            Sequence = proteinXmlEntry.Sequence;
            SourceDatabase = PeffEntrySourceDatabase.Uniprot; // Default to Uniprot, can be changed later
        }
        public ProteinPeffEntry(Protein protein, PeffEntrySourceDatabase source = PeffEntrySourceDatabase.Uniprot)
        {
            Accession = protein.Accession;
            //Description = protein.Description;
            //GeneName = protein.GeneName;
            Organism = protein.Organism;
            Sequence = protein.BaseSequence;
            SourceDatabase = source;
        }
        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append($">{PeffEntrySourceDatabaseAbbreviation(SourceDatabase)}:{Accession}");
            sb.Append($" \\DbUniqueId={Accession}");
            sb.Append($" \\OX={Organism}");
            sb.AppendLine($" \\Length= {Sequence.Length}");
            int currentPosition = 0;
            int lineLength = 60; // Length of each line in the wrapped string
            StringBuilder wrappedString = new StringBuilder();
            while (currentPosition < Sequence.Length)
            {
                int remainingLength = Sequence.Length - currentPosition;
                int lengthToTake = Math.Min(lineLength, remainingLength);

                wrappedString.Append(Sequence.Substring(currentPosition, lengthToTake));

                currentPosition += lengthToTake;

                if (currentPosition < Sequence.Length) // Add a newline if there's more text
                {
                    wrappedString.Append(Environment.NewLine);
                }
            }
            sb.AppendLine(wrappedString.ToString());
            return sb.ToString();
        }
        public enum PeffEntrySourceDatabase
        {
            Uniprot,
            Ensembl,
            RefSeq,
            NCBI,
            Custom
        }
        public static string PeffEntrySourceDatabaseAbbreviation(PeffEntrySourceDatabase source) => source switch
        {
            PeffEntrySourceDatabase.Uniprot => "up",
            PeffEntrySourceDatabase.Ensembl => "es",
            PeffEntrySourceDatabase.RefSeq => "",
            PeffEntrySourceDatabase.NCBI => "nc",
            _ => ""
        };
    }
}
