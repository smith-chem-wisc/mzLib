﻿using System;
using System.Collections.Generic;
using System.Linq;

namespace Proteomics
{
    public class Protein
    {
        #region Public Constructors

        public Protein(string sequence, string accession, string organism = null, List<Tuple<string, string>> gene_names = null, IDictionary<int, List<Modification>> oneBasedModifications = null, List<ProteolysisProduct> proteolysisProducts = null, string name = null, string full_name = null, bool isDecoy = false, bool isContaminant = false, List<DatabaseReference> databaseReferences = null, List<SequenceVariation> sequenceVariations = null, List<DisulfideBond> disulfideBonds = null, string databaseFilePath = null)
        {
            // Mandatory
            BaseSequence = sequence;
            Accession = accession;

            Name = name;
            Organism = organism;
            FullName = full_name;
            IsDecoy = isDecoy;
            IsContaminant = isContaminant;
            DatabaseFilePath = databaseFilePath;

            GeneNames = gene_names ?? new List<Tuple<string, string>>();
            ProteolysisProducts = proteolysisProducts ?? new List<ProteolysisProduct>();
            SequenceVariations = sequenceVariations ?? new List<SequenceVariation>();
            OneBasedPossibleLocalizedModifications = oneBasedModifications ?? new Dictionary<int, List<Modification>>();
            DatabaseReferences = databaseReferences ?? new List<DatabaseReference>();
            DisulfideBonds = disulfideBonds ?? new List<DisulfideBond>();
        }

        #endregion Public Constructors

        #region Public Properties

        public IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; }

        /// <summary>
        /// The list of gene names consists of tuples, where Item1 is the type of gene name, and Item2 is the name. There may be many genes and names of a certain type produced when reading an XML protein database.
        /// </summary>
        public IEnumerable<Tuple<string, string>> GeneNames { get; }

        public string Accession { get; }
        public string BaseSequence { get; }
        public string Organism { get; }
        public bool IsDecoy { get; }
        public IEnumerable<SequenceVariation> SequenceVariations { get; }
        public IEnumerable<DisulfideBond> DisulfideBonds { get; }
        public IEnumerable<ProteolysisProduct> ProteolysisProducts { get; }
        public IEnumerable<DatabaseReference> DatabaseReferences { get; }
        public string DatabaseFilePath { get; }

        public int Length
        {
            get
            {
                return BaseSequence.Length;
            }
        }

        public string FullDescription
        {
            get
            {
                return Accession + "|" + Name + "|" + FullName;
            }
        }

        public string Name { get; }

        public string FullName { get; }

        public bool IsContaminant { get; }

        #endregion Public Properties

        #region Public Indexers

        public char this[int zeroBasedIndex]
        {
            get
            {
                return BaseSequence[zeroBasedIndex];
            }
        }

        #endregion Public Indexers

        #region Public Methods

        /// <summary>
        /// Formats a string for a UniProt fasta header. See https://www.uniprot.org/help/fasta-headers.
        /// Note that the db field isn't very applicable here, so mz is placed in to denote written by mzLib.
        /// </summary>
        /// <returns></returns>
        public string GetUniProtFastaHeader()
        {
            var n = GeneNames.FirstOrDefault();
            string geneName = n == null ? "" : n.Item2;
            return String.Format("mz|{0}|{1} {2} OS={3} GN={4}", Accession, Name, FullName, Organism, geneName);
        }

        /// <summary>
        /// Formats a string for an ensembl header
        /// </summary>
        /// <returns></returns>
        public string GetEnsemblFastaHeader()
        {
            return String.Format("{0} {1}", Accession, FullName);
        }

        #endregion Public Methods
    }
}