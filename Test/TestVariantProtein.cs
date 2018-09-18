using NUnit.Framework;
using Proteomics;
using System.Collections.Generic;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public class VariantProteinTests
    {
        [Test]
        public void MakeVariantFasta()
        {
            List<Protein> proteins = ProteinDbLoader.LoadProteinXML(@"E:\ProjectsActive\JurkatProteogenomics\180413\combined_1-trimmed-pair1Aligned.sortedByCoord.outProcessed.out.fixedQuals.split.concat.sorted.snpEffAnnotated.protein.xml",
                true, DecoyType.None, null, false, null, out var un);
            List<ProteinWithAppliedVariants> variantProteins = proteins.SelectMany(p => p.GetVariantProteins()).ToList();
            List<ProteinWithAppliedVariants> variantProteins2 = variantProteins.Select(p =>
                new ProteinWithAppliedVariants(
                    p.BaseSequence,
                    new Protein(p.Protein.BaseSequence, p.Protein.Accession, p.Protein.Organism, p.Protein.GeneNames.ToList(), p.Protein.OneBasedPossibleLocalizedModifications, p.Protein.ProteolysisProducts.ToList(), p.Protein.Name + string.Join(",", p.AppliedSequenceVariations.Select(d => d.Description)), p.Protein.FullName + string.Join(",", p.AppliedSequenceVariations.Select(d => d.Description)), p.IsDecoy, p.IsContaminant, p.DatabaseReferences.ToList(), p.SequenceVariations.ToList(), p.DisulfideBonds.ToList(), p.DatabaseFilePath),
                    p.AppliedSequenceVariations, p.Individual)
                ).ToList();
            ProteinDbWriter.WriteFastaDatabase(variantProteins2.OfType<Protein>().ToList(), @"E:\ProjectsActive\Spritz\mmTesting\variantproteins.fasta", "|");
        }
    }
}