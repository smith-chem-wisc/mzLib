using System.Collections.Generic;
using NUnit.Framework;
using Omics.Digestion;
using Assert = NUnit.Framework.Legacy.ClassicAssert;

namespace Test.Omics.Modifications
{
    /// <summary>
    /// Pins DigestionMotif's Schechter-Berger subsite accessors against motifs that ship, including
    /// motifs that straddle the bond (collagenase "GPX|GPX", StcE-trypsin "TX|T", colicin_E5 "G|U").
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class DigestionMotifSubsiteTests
    {
        private static DigestionMotif Motif(string motifString)
        {
            List<DigestionMotif> parsed = DigestionMotif.ParseDigestionMotifsFromString(motifString);
            Assert.AreEqual(1, parsed.Count, $"test setup: '{motifString}' should parse to one motif");
            return parsed[0];
        }

        [Test]
        public static void Subsites_Trypsin_NamesP1Only()
        {
            DigestionMotif trypsin = Motif("K|");
            Assert.AreEqual('K', trypsin.NonPrimeSubsite(1));
            Assert.AreEqual('\0', trypsin.NonPrimeSubsite(2));
            Assert.AreEqual('\0', trypsin.PrimeSubsite(1), "trypsin's motif says nothing about P1'");
        }

        [Test]
        public static void Subsites_AspN_NamesP1PrimeOnly()
        {
            DigestionMotif aspN = Motif("|D");
            Assert.AreEqual('\0', aspN.NonPrimeSubsite(1));
            Assert.AreEqual('D', aspN.PrimeSubsite(1));
            Assert.AreEqual('\0', aspN.PrimeSubsite(2));
        }

        [Test]
        public static void Subsites_StcETrypsin_NamesBothSides_WithTheGlycosylatableResidueAtP2()
        {
            // "TX|T" -- the literature places StcE's glycan requirement at P2, and P2 is exactly where
            // this motif names a Thr.
            DigestionMotif stcE = Motif("TX|T");
            Assert.AreEqual('X', stcE.NonPrimeSubsite(1));
            Assert.AreEqual('T', stcE.NonPrimeSubsite(2));
            Assert.AreEqual('\0', stcE.NonPrimeSubsite(3));
            Assert.AreEqual('T', stcE.PrimeSubsite(1));
            Assert.AreEqual('\0', stcE.PrimeSubsite(2));
        }

        [Test]
        public static void Subsites_Collagenase_WalksBothDirectionsAcrossTheBond()
        {
            DigestionMotif collagenase = Motif("GPX|GPX");
            Assert.AreEqual('X', collagenase.NonPrimeSubsite(1));
            Assert.AreEqual('P', collagenase.NonPrimeSubsite(2));
            Assert.AreEqual('G', collagenase.NonPrimeSubsite(3));
            Assert.AreEqual('\0', collagenase.NonPrimeSubsite(4));
            Assert.AreEqual('G', collagenase.PrimeSubsite(1));
            Assert.AreEqual('P', collagenase.PrimeSubsite(2));
            Assert.AreEqual('X', collagenase.PrimeSubsite(3));
            Assert.AreEqual('\0', collagenase.PrimeSubsite(4));
        }

        [Test]
        public static void Subsites_ColicinE5_OneResidueEitherSide()
        {
            DigestionMotif colicin = Motif("G|U");
            Assert.AreEqual('G', colicin.NonPrimeSubsite(1));
            Assert.AreEqual('U', colicin.PrimeSubsite(1));
        }

        [Test]
        [TestCase(0)]
        [TestCase(-1)]
        public static void Subsites_PositionBelowOne_IsNotASubsite(int position)
        {
            // P0 does not exist; without the guard NonPrimeSubsite(0) would return the P1' residue.
            DigestionMotif stcE = Motif("TX|T");
            Assert.AreEqual('\0', stcE.NonPrimeSubsite(position));
            Assert.AreEqual('\0', stcE.PrimeSubsite(position));
        }

        [Test]
        public static void Subsites_NoRecognitionSequence_NamesNoSubsite()
        {
            // The public constructor accepts a null recognition sequence.
            DigestionMotif motif = new DigestionMotif(null, null, 1, null);
            Assert.AreEqual('\0', motif.NonPrimeSubsite(1));
            Assert.AreEqual('\0', motif.PrimeSubsite(1));
        }
    }
}
