using FlashLFQ;
using NUnit.Framework;
using Omics;
using Omics.BioPolymerGroup;
using Omics.Modifications;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Transcriptomics;

namespace Test.Omics
{
    [ExcludeFromCodeCoverage]
    internal class BioPolymerGroupTests
    {
        [Test]
        public static void TestProteinGroupsAccessionOutputOrder()
        {
            var biopolymers = new HashSet<IBioPolymer>();
            var biopolymerWithSetMods = new HashSet <IBioPolymerWithSetMods>();
            var uniqueBioPolymerWithSetMods = new HashSet< IBioPolymerWithSetMods >();
            
            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();

            // make protein B
            biopolymers.Add(new Protein("-----F----*", "B", null, gn, new Dictionary<int, List<Modification>>(), isDecoy: true));

            // make protein A
            biopolymers.Add(new Protein("-----F----**", "A", null, gn, new Dictionary<int, List<Modification>>(), isDecoy: true));

            // add protein B and A to the protein group
            BioPolymerGroup testGroup = new BioPolymerGroup(biopolymers, biopolymerWithSetMods, uniqueBioPolymerWithSetMods);

            // test order is AB and not BA
            Assert.That(testGroup.BioPolymerGroupName.Equals("A|B"));
            Assert.That(testGroup.BioPolymers.First().Accession.Equals("B"));
        }
        [Test]
        public static void TestRNAGroupsAccessionOutputOrder()
        {
            var biopolymers = new HashSet<IBioPolymer>();
            var biopolymerWithSetMods = new HashSet<IBioPolymerWithSetMods>();
            var uniqueBioPolymerWithSetMods = new HashSet<IBioPolymerWithSetMods>();

            List<Tuple<string, string>> gn = new List<Tuple<string, string>>();

            // make protein B
            biopolymers.Add(new RNA("AAAAACAAAAU", "B", isDecoy: true));

            // make protein A
            biopolymers.Add(new RNA("AAAAACAAAAUU", "A", isDecoy: true));

            // add protein B and A to the protein group
            BioPolymerGroup testGroup = new BioPolymerGroup(biopolymers, biopolymerWithSetMods, uniqueBioPolymerWithSetMods);

            // test order is AB and not BA
            Assert.That(testGroup.BioPolymerGroupName.Equals("A|B"));
            Assert.That(testGroup.BioPolymers.First().Accession.Equals("B"));
        }
    }
}
