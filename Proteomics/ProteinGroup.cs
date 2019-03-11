using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Proteomics
{
    public class ProteinGroup : IEquatable<ProteinGroup>
    {
        public ProteinList ProteinList { get; private set; }
        public ParsimonySequenceList PepWithSetModsList { get; private set; }

        public ProteinGroup(ProteinList proteinList = null, ParsimonySequenceList pwsmList = null)

        {
            this.ProteinList = new ProteinList();
            this.PepWithSetModsList = new ParsimonySequenceList();

            if (proteinList != null)
            {
                this.ProteinList.IntegrateProteinLists(proteinList);
            }
            if (pwsmList != null)
            {
                this.PepWithSetModsList.AddParsimonySequences(pwsmList);
            }
        }

        public void AddToThisProteinGroup(ProteinGroup group)
        {
            this.ProteinList.IntegrateProteinLists(group.ProteinList);
            this.PepWithSetModsList.AddParsimonySequences(group.PepWithSetModsList);
        }

        public bool Equals(ProteinGroup that)
        {
            return that != null && that.ProteinList == this.ProteinList && that.PepWithSetModsList == this.PepWithSetModsList;
        }

        public override bool Equals(object that)
        {
            return this.Equals(that as ProteinGroup);
        }

        public override int GetHashCode()
        {
            this.PepWithSetModsList.ParsimonySequences.Sort();

            var aHashCode = this.ProteinList.Proteins[0].GetHashCode();

            for (int i = 1; i < this.ProteinList.Proteins.Count; i++)
            {
                aHashCode = aHashCode ^ this.ProteinList.Proteins[i].GetHashCode();
            }

            aHashCode = aHashCode ^ this.PepWithSetModsList.ParsimonySequences[0].GetHashCode();
            for (int i = 1; i < this.PepWithSetModsList.ParsimonySequences.Count; i++)
            {
                aHashCode = aHashCode ^ this.PepWithSetModsList.ParsimonySequences[i].GetHashCode();
            }
            return aHashCode;
        }

        public override string ToString()
        {
            return string.Join(", ", ProteinList.AccessionNumbers);
        }
    }

    public class ProteinList
    {
        public List<Protein> Proteins { get; private set; }
        public List<string> AccessionNumbers { get; }

        public ProteinList(List<Protein> proteinList = null)
        {
            this.Proteins = new List<Protein>();
            this.AccessionNumbers = new List<string>();

            if (proteinList != null)
            {
                this.Proteins.AddRange(proteinList.Distinct());
                this.AccessionNumbers.AddRange(proteinList.Select(p => p.Accession).Distinct());
                Proteins.Sort((x, y) => x.Accession.CompareTo(y.Accession));
                AccessionNumbers.Sort((x, y) => x.CompareTo(y));
            }
        }

        public void AddProtein(Protein newProtein)
        {
            if (!this.Proteins.Contains(newProtein))
            {
                this.Proteins.Add(newProtein);
                this.AccessionNumbers.Add(newProtein.Accession);
            }

            //TODO: make this more efficient. this is re-sorting the entire list every time a new protein gets added
            Proteins.Sort((x, y) => x.Accession.CompareTo(y.Accession));
            AccessionNumbers.Sort((x, y) => x.CompareTo(y));
        }

        public void IntegrateProteinLists(ProteinList newProteinsToAdd)
        {
            var accessionsToAdd = newProteinsToAdd.AccessionNumbers.Except(this.AccessionNumbers);

            foreach (Protein p in newProteinsToAdd.Proteins.Distinct())
            {
                if (accessionsToAdd.Contains(p.Accession))
                {
                    this.Proteins.Add(p);
                    this.AccessionNumbers.Add(p.Accession);
                }
            }
        }

        public override bool Equals(object other)
        {
            ProteinList otherList = (ProteinList)other;

            if (otherList == null)
            {
                return false;
            }

            return otherList.Proteins.SequenceEqual(this.Proteins);
        }

        public override int GetHashCode()
        {
            if (this.Proteins == null || !this.Proteins.Any())
            {
                // maybe change this?
                return 0;
                //return base.GetHashCode();
            }

            return this.Proteins[0].Accession.GetHashCode();
        }

        public override string ToString()
        {
            return String.Join(", ", this.AccessionNumbers).ToString();
        }
    }

    public class ParsimonySequenceList
    {
        //TODO: make sure these lists are sorted, same style as proteinList
        public List<ParsimonySequence> ParsimonySequences { get; set; }

        public List<string> FullSequenceList { get; }

        public ParsimonySequenceList(List<ParsimonySequence> pwsmList = null)
        {
            this.ParsimonySequences = new List<ParsimonySequence>();
            this.FullSequenceList = new List<string>();

            if (pwsmList != null)
            {
                this.ParsimonySequences.AddRange(pwsmList);
                this.FullSequenceList.AddRange(pwsmList.Select(s => s.Sequence).ToList());
            }
        }

        public void AddParsimonySequences(ParsimonySequenceList newParsimonySequenceList)
        {
            if (this.ParsimonySequences.Count == 0 && newParsimonySequenceList != null)
            {
                foreach (ParsimonySequence p in newParsimonySequenceList.ParsimonySequences.Distinct())
                {
                    this.ParsimonySequences.Add(p);
                    this.FullSequenceList.Add(p.Sequence);
                }
            }
            else
            {
                var parsimonySequencesToAdd =  newParsimonySequenceList.ParsimonySequences.Except(this.ParsimonySequences).ToList();
                foreach (ParsimonySequence p in parsimonySequencesToAdd.Distinct())
                {
                    this.ParsimonySequences.Add(p);
                    this.FullSequenceList.Add(p.Sequence);
                }
            }
        }

        public override bool Equals(object other)
        {
            ParsimonySequenceList otherList = (ParsimonySequenceList)other;

            if (otherList == null)
            {
                return false;
            }

            return otherList.ParsimonySequences.SequenceEqual(this.ParsimonySequences);
        }

        public override int GetHashCode()
        {
            if (this.ParsimonySequences == null || !this.ParsimonySequences.Any())
            {
                // maybe change this?
                return 0;
                //return base.GetHashCode();
            }

            return this.ParsimonySequences[0].Sequence.GetHashCode();
        }

        public override string ToString()
        {
            return string.Join(", ", FullSequenceList);
        }
    }
}