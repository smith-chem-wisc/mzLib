using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using MassSpectrometry;
using Proteomics.ProteolyticDigestion;
using Proteomics;
using NUnit.Framework;

namespace Test
{
    [TestFixture]
    public static class TestIsotopicEnvelope
    {
        [Test]
        public static void TestIsotopicDistributionToEnvelope()
        {
            Protein myProtein = new Protein("PEPTIDEK", "accession");

            DigestionParams digest = new DigestionParams(protease: "trypsin", maxMissedCleavages: 0,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
           
            PeptideWithSetModifications pwsm = myProtein.Digest(digest, new List<Modification>(),
                new List<Modification>()).First();
           
            IsotopicDistribution distribution = IsotopicDistribution.GetDistribution(pwsm.FullChemicalFormula, 0.125, 1e-8);
            IsotopicEnvelope envelope = new(distribution, charge: 2);

            double distributionMostAbundant = distribution.MostAbundantMass;
            double envelopeMostAbundant = envelope.MostAbundantObservedIsotopicMass;
            Assert.AreEqual(envelopeMostAbundant, envelopeMostAbundant);
        }

    }
}
