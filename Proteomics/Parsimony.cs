using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Proteomics
{
    public class Parsimony
    {
        public List<ProteinGroup> ProteinGroups { get; private set; }

        public Parsimony(List<ProteinGroup> proteinGroups = null)
        {
            this.ProteinGroups = new List<ProteinGroup>();
            if (proteinGroups != null)
            {
                this.ProteinGroups.AddRange(proteinGroups);
            }
        }

        public void AssignProteingroups(List<PeptideWithSetModifications> observedPeptides, List<Protein> theDatabase)
        {
            //This dictionary has a list of all proteins that could be the origin of a theoretical peptide and is formed from the complete theoretical database
            Dictionary<ParsimonySequence, ProteinList> TheoreticalParsimonySequenceToProteinLookupTable = GetTheoreticalPeptideAccessionDictionary(theDatabase);

            //This dictioary takes all the experimentally observed peptides and determines which proteins they could have come from based on the dictionary above.
            //this is a comprehensive/exhaustive set
            List<Tuple<ParsimonySequenceList, ProteinList>> ExperimentalParsimonySequenceToProteinLookupTable = GetExperimentalPeptideAccessionDictionary(TheoreticalParsimonySequenceToProteinLookupTable, observedPeptides);

            ////////////////////////////////////////////////////////////////////////////////////////////////////
            //get all the unique protein lists
            var possibleProteins = ExperimentalParsimonySequenceToProteinLookupTable.Select(p => p.Item2).Distinct().ToList();

            //foreach protein in the proteinlist
            foreach (ProteinList proteinCandidate in possibleProteins)
            {
                //get indicies all the peptides that match to that protein
                var indexes = new List<int>();
                for (int i = 0; i < ExperimentalParsimonySequenceToProteinLookupTable.Count; i++)
                {
                    if (ExperimentalParsimonySequenceToProteinLookupTable[i].Item2.Equals(proteinCandidate))
                    {
                        indexes.Add(i);
                    }
                }

                if (indexes.Count > 1)
                {
                    if (ExperimentalParsimonySequenceToProteinLookupTable.Where(c => c.Item1.Equals(ExperimentalParsimonySequenceToProteinLookupTable[indexes[0]].Item1)).Count() == 1)
                    {
                        for (int i = indexes.Count - 1; i > 0; i--)
                        {
                            if (ExperimentalParsimonySequenceToProteinLookupTable.Where(c => c.Item1.Equals(ExperimentalParsimonySequenceToProteinLookupTable[indexes[i]].Item1)).Count() == 1)
                            {
                                ExperimentalParsimonySequenceToProteinLookupTable[indexes[0]].Item1.FullSequenceList.AddRange(ExperimentalParsimonySequenceToProteinLookupTable[indexes[i]].Item1.FullSequenceList);
                                ExperimentalParsimonySequenceToProteinLookupTable[indexes[0]].Item1.ParsimonySequences.AddRange(ExperimentalParsimonySequenceToProteinLookupTable[indexes[i]].Item1.ParsimonySequences);
                                ExperimentalParsimonySequenceToProteinLookupTable.RemoveAt(indexes[i]);
                            }
                        }
                    }
                }
            }

            ////////////////////////////////////////////////////////////////////////////////////////

            Dictionary<ProteinList, List<ParsimonySequenceList>> d = new Dictionary<ProteinList, List<ParsimonySequenceList>>();
            foreach (ProteinList proteinCandidate in possibleProteins)
            {
                //get indicies all the peptides that match to that protein
                List<ParsimonySequenceList> pslList = new List<ParsimonySequenceList>();
                for (int i = 0; i < ExperimentalParsimonySequenceToProteinLookupTable.Count; i++)
                {
                    if (ExperimentalParsimonySequenceToProteinLookupTable[i].Item2.Equals(proteinCandidate))
                    {
                        if (d.ContainsKey(proteinCandidate))
                        {
                            d[proteinCandidate].Add(ExperimentalParsimonySequenceToProteinLookupTable[i].Item1);
                        }
                        else
                        {
                            d.Add(proteinCandidate, new List<ParsimonySequenceList>() { ExperimentalParsimonySequenceToProteinLookupTable[i].Item1 });
                        }
                    }
                }
            }

            List<List<ProteinList>> llp = new List<List<ProteinList>>();
            List<ProteinList> ignore = new List<ProteinList>();
            foreach (KeyValuePair<ProteinList, List<ParsimonySequenceList>> kvp in d)
            {
                if (!ignore.Contains(kvp.Key))
                {
                    List<ProteinList> lp = new List<ProteinList>() { kvp.Key };

                    foreach (KeyValuePair<ProteinList, List<ParsimonySequenceList>> kvp2 in d)
                    {
                        if (kvp2.Key != kvp.Key)
                        {
                            if (kvp2.Value.Except(kvp.Value).Count() == 0 && kvp.Value.Except(kvp2.Value).Count() == 0)
                            {
                                lp.Add(kvp2.Key);
                            }
                        }
                    }

                    if (lp.Count > 1)
                    {
                        llp.Add(lp);
                        foreach (ProteinList item in lp)
                        {
                            ignore.Add(item);
                        }
                    }
                }
            }



            List<Tuple<ParsimonySequenceList, ProteinList>> e = new List<Tuple<ParsimonySequenceList, ProteinList>>();
            foreach (var item in llp)
            {
                ProteinList myPl = new ProteinList();
                foreach (ProteinList p in item)
                {
                    myPl.AccessionNumbers.AddRange(p.AccessionNumbers);
                    myPl.Proteins.AddRange(p.Proteins);
                }


                for (int i = ExperimentalParsimonySequenceToProteinLookupTable.Count - 1; i >= 0; i--)
                {
                    if (item.Contains(ExperimentalParsimonySequenceToProteinLookupTable[i].Item2))
                    {
                        Tuple<ParsimonySequenceList, ProteinList> tp = new Tuple<ParsimonySequenceList, ProteinList>(ExperimentalParsimonySequenceToProteinLookupTable[i].Item1, myPl);
                        ExperimentalParsimonySequenceToProteinLookupTable.RemoveAt(i);
                        e.Add(tp);
                    }
                }

            }

            ExperimentalParsimonySequenceToProteinLookupTable.AddRange(e);
        }

        private List<Tuple<ParsimonySequenceList, ProteinList>> GetExperimentalPeptideAccessionDictionary(Dictionary<ParsimonySequence, ProteinList> theoreticalParsimonySequenceToProteinLookupTable, List<PeptideWithSetModifications> observedPeptides)
        {
            List<Tuple<ParsimonySequenceList, ProteinList>> d = new List<Tuple<ParsimonySequenceList, ProteinList>>();
            foreach (PeptideWithSetModifications peptide in observedPeptides)
            {
                //TODO: set the boolean based on user settings for treat modified peptides as unique
                ParsimonySequence ps = new ParsimonySequence(peptide, false);
                ParsimonySequenceList pl = new ParsimonySequenceList(new List<ParsimonySequence>() { ps });

                foreach (Protein p in theoreticalParsimonySequenceToProteinLookupTable[ps].Proteins)
                {
                    ProteinList newPl = new ProteinList(new List<Protein>() { p });
                    Tuple<ParsimonySequenceList, ProteinList> t = new Tuple<ParsimonySequenceList, ProteinList>(pl, newPl);

                    if (!d.Contains(t))
                    {
                        d.Add(t);
                    }
                }
            }
            return d;
        }

        private List<ProteinGroup> CombineGroupsWithSameProteinList(List<ProteinGroup> myProteinGroups)
        {
            List<ProteinList> allTheDifferentProteinLists = new List<ProteinList>();
            foreach (ProteinGroup pg in myProteinGroups)
            {
                if (!allTheDifferentProteinLists.Contains(pg.ProteinList))
                {
                    allTheDifferentProteinLists.Add(pg.ProteinList);
                }
            }

            List<ProteinGroup> returnProteinGroupList = new List<ProteinGroup>();

            foreach (ProteinList proteinList in allTheDifferentProteinLists)
            {
                ProteinGroup newProteinGroup = new ProteinGroup();
                foreach (ProteinGroup pg in myProteinGroups.Where(pl => pl.ProteinList.Equals(proteinList)).ToList())
                {
                    newProteinGroup.AddToThisProteinGroup(pg);
                }
                returnProteinGroupList.Add(newProteinGroup);
            }

            return returnProteinGroupList;
        }

        private List<ProteinGroup> CombineProteinGroupsWherePeptidesMatchToTheSameProteinList(List<ProteinGroup> myProteinGroups)
        {
            List<ParsimonySequenceList> allTheDifferentPeptideLists = new List<ParsimonySequenceList>();
            foreach (ProteinGroup pg in myProteinGroups)
            {
                if (!allTheDifferentPeptideLists.Contains(pg.PepWithSetModsList))
                {
                    allTheDifferentPeptideLists.Add(pg.PepWithSetModsList);
                }
            }

            List<ProteinGroup> returnProteinGroupList = new List<ProteinGroup>();

            foreach (ParsimonySequenceList peptideList in allTheDifferentPeptideLists)
            {
                List<ProteinGroup> tempProteinGroupListOne = new List<ProteinGroup>();
                tempProteinGroupListOne.AddRange(myProteinGroups.Where(p => p.PepWithSetModsList.Equals(peptideList)));
            }

            foreach (ParsimonySequenceList peptideList in allTheDifferentPeptideLists)
            {
                ProteinGroup newProteinGroup = new ProteinGroup();

                foreach (ProteinGroup pg in myProteinGroups.Where(p => p.PepWithSetModsList.Equals(peptideList)).ToList())
                {
                    newProteinGroup.AddToThisProteinGroup(pg);
                }
                returnProteinGroupList.Add(newProteinGroup);
            }

            return returnProteinGroupList;
        }

        private void ShrinkTheSet()
        {
            int beforeCount = ProteinGroups.Count;
            int afterCount = 0;
            do
            {
                beforeCount = ProteinGroups.Count;

                this.ProteinGroups = CombineProteinGroupsWherePeptidesMatchToTheSameProteinList(this.ProteinGroups);
                this.ProteinGroups = CombineGroupsWithSameProteinList(this.ProteinGroups);

                afterCount = ProteinGroups.Count;
            } while (beforeCount > afterCount);
        }

        private Dictionary<ParsimonySequence, ProteinList> GetTheoreticalPeptideAccessionDictionary(List<Protein> theDatabase)
        {
            Dictionary<ParsimonySequence, ProteinList> peptideAccessionDictionary = new Dictionary<ParsimonySequence, ProteinList>();

            foreach (Protein protein in theDatabase)
            {
                //TODO: the digestion params should match that used for the sample.
                var digestedProtein = protein.Digest(new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1), new List<Modification>(), new List<Modification>());

                foreach (PeptideWithSetModifications peptide in digestedProtein)
                {
                    //TODO: the boolean needs to match the user selection (treat modified peptides as unique)
                    ParsimonySequence ps = new ParsimonySequence(peptide, false);

                    if (peptideAccessionDictionary.ContainsKey(ps))
                    {
                        if (!peptideAccessionDictionary[ps].Proteins.Contains(protein))
                        {
                            peptideAccessionDictionary[ps].AddProtein(protein);
                        }
                    }
                    else
                    {
                        peptideAccessionDictionary.Add(ps, new ProteinList(new List<Protein>() { protein }));
                    }
                }
            }

            return peptideAccessionDictionary;
        }
    }
}