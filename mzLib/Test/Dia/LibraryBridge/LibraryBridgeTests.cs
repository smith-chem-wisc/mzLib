// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;
using MassSpectrometry.Dia;
using MzLibUtil;
using NUnit.Framework;

namespace Test.Dia.LibraryBridge
{
    [TestFixture]
    public class LibraryBridgeTests
    {
        private DiaScanIndex _scanIndex;

        /// <summary>
        /// Builds a DiaScanIndex with 3 DIA windows via DiaScanIndexBuilder.
        ///   Window 0: center 412.5, width 25 → [400, 425]
        ///   Window 1: center 437.5, width 25 → [425, 450]
        ///   Window 2: center 475.0, width 50 → [450, 500]
        /// 2 cycles × 3 windows = 6 MS2 scans, RT 10–13 min.
        /// </summary>
        [SetUp]
        public void Setup()
        {
            var scans = new List<MsDataScan>();
            int scanNumber = 1;

            var windows = new[]
            {
                (center: 412.5, width: 25.0),
                (center: 437.5, width: 25.0),
                (center: 475.0, width: 50.0),
            };

            double[] cycleStartRts = { 10.0, 12.0 };

            foreach (var cycleRt in cycleStartRts)
            {
                for (int w = 0; w < windows.Length; w++)
                {
                    double rt = cycleRt + w * 0.5;
                    double[] mzs = { 100.1, 200.2, 300.3, 400.4, 500.5 };
                    double[] intensities = { 1000, 2000, 3000, 4000, 5000 };

                    var spectrum = new MzSpectrum(mzs, intensities, false);
                    var scan = new MsDataScan(
                        massSpectrum: spectrum,
                        oneBasedScanNumber: scanNumber++,
                        msnOrder: 2,
                        isCentroid: true,
                        polarity: Polarity.Positive,
                        retentionTime: rt,
                        scanWindowRange: new MzRange(100, 600),
                        scanFilter: "",
                        mzAnalyzer: MZAnalyzerType.Orbitrap,
                        totalIonCurrent: intensities.Sum(),
                        injectionTime: 50,
                        noiseData: null,
                        nativeId: $"scan={scanNumber - 1}",
                        isolationMZ: windows[w].center,
                        isolationWidth: windows[w].width,
                        dissociationType: DissociationType.HCD
                    );
                    scans.Add(scan);
                }
            }

            _scanIndex = DiaScanIndexBuilder.Build(scans.ToArray());
        }

        private static LibraryPrecursorInput MakePrecursor(
            string sequence, double precursorMz, int charge,
            float[] fragmentMzs, float[] fragmentIntensities,
            double? rt = null, bool isDecoy = false)
        {
            return new LibraryPrecursorInput(
                sequence, precursorMz, charge, rt, isDecoy,
                fragmentMzs, fragmentIntensities);
        }

        // ═══════════════════════════════════════════════════════════════════
        // FindWindowForPrecursorMz
        // ═══════════════════════════════════════════════════════════════════

        [Test]
        public void FindWindow_PrecursorInFirstWindow()
        {
            Assert.That(_scanIndex.FindWindowForPrecursorMz(410f), Is.EqualTo(0));
        }

        [Test]
        public void FindWindow_PrecursorInMiddleWindow()
        {
            Assert.That(_scanIndex.FindWindowForPrecursorMz(435f), Is.EqualTo(1));
        }

        [Test]
        public void FindWindow_PrecursorInLastWindow()
        {
            Assert.That(_scanIndex.FindWindowForPrecursorMz(475f), Is.EqualTo(2));
        }

        [Test]
        public void FindWindow_BelowAllWindows_ReturnsNegativeOne()
        {
            Assert.That(_scanIndex.FindWindowForPrecursorMz(350f), Is.EqualTo(-1));
        }

        [Test]
        public void FindWindow_AboveAllWindows_ReturnsNegativeOne()
        {
            Assert.That(_scanIndex.FindWindowForPrecursorMz(550f), Is.EqualTo(-1));
        }

        [Test]
        public void FindWindow_DoubleOverload()
        {
            Assert.That(_scanIndex.FindWindowForPrecursorMz(435.0), Is.EqualTo(1));
            Assert.That(_scanIndex.FindWindowForPrecursorMz(350.0), Is.EqualTo(-1));
        }

        [Test]
        public void FindWindow_AtBoundary_ReturnsValidWindow()
        {
            int result = _scanIndex.FindWindowForPrecursorMz(425f);
            Assert.That(result, Is.GreaterThanOrEqualTo(0));
            Assert.That(result, Is.LessThan(_scanIndex.WindowCount));
        }

        [Test]
        public void GetGlobalRtMin_ReturnsMinRt()
        {
            Assert.That(_scanIndex.GetGlobalRtMin(), Is.EqualTo(10.0f).Within(0.01f));
        }

        [Test]
        public void GetGlobalRtMax_ReturnsMaxRt()
        {
            Assert.That(_scanIndex.GetGlobalRtMax(), Is.EqualTo(13.0f).Within(0.01f));
        }

        // ═══════════════════════════════════════════════════════════════════
        // DiaLibraryQueryGenerator.Generate()
        // ═══════════════════════════════════════════════════════════════════

        [Test]
        public void Generate_SinglePrecursor_CorrectQueryCount()
        {
            var p = new DiaSearchParameters { PpmTolerance = 20f, RtToleranceMinutes = 2f, MinFragmentsRequired = 1 };
            var precursors = new List<LibraryPrecursorInput>
            {
                MakePrecursor("PEPTIDE", 412.5, 2, new float[] { 200.1f, 300.2f, 400.3f }, new float[] { 100f, 80f, 60f }, rt: 11.0)
            };
            var result = DiaLibraryQueryGenerator.Generate(precursors, _scanIndex, p);
            Assert.That(result.Queries.Length, Is.EqualTo(3));
            Assert.That(result.PrecursorGroups.Length, Is.EqualTo(1));
            Assert.That(result.SkippedNoWindow, Is.EqualTo(0));
        }

        [Test]
        public void Generate_QueryFieldsCorrect()
        {
            var p = new DiaSearchParameters { PpmTolerance = 20f, RtToleranceMinutes = 2f, MinFragmentsRequired = 1 };
            var precursors = new List<LibraryPrecursorInput>
            {
                MakePrecursor("PEPTIDE", 412.5, 2, new float[] { 200.1f }, new float[] { 100f }, rt: 11.0)
            };
            var result = DiaLibraryQueryGenerator.Generate(precursors, _scanIndex, p);
            var q = result.Queries[0];
            Assert.That(q.TargetMz, Is.EqualTo(200.1f).Within(0.01f));
            Assert.That(q.TolerancePpm, Is.EqualTo(20f));
            Assert.That(q.RtMin, Is.EqualTo(9.0f).Within(0.01f));
            Assert.That(q.RtMax, Is.EqualTo(13.0f).Within(0.01f));
            Assert.That(q.WindowId, Is.EqualTo(0));
        }

        [Test]
        public void Generate_MultiplePrecursors_QueriesContiguous()
        {
            var p = new DiaSearchParameters { PpmTolerance = 20f, RtToleranceMinutes = 2f, MinFragmentsRequired = 1 };
            var precursors = new List<LibraryPrecursorInput>
            {
                MakePrecursor("A", 412.5, 2, new float[] { 200f, 300f }, new float[] { 100f, 80f }, rt: 11.0),
                MakePrecursor("B", 437.5, 2, new float[] { 250f, 350f, 450f }, new float[] { 90f, 70f, 50f }, rt: 12.0)
            };
            var result = DiaLibraryQueryGenerator.Generate(precursors, _scanIndex, p);
            Assert.That(result.Queries.Length, Is.EqualTo(5));
            Assert.That(result.PrecursorGroups.Length, Is.EqualTo(2));
            Assert.That(result.PrecursorGroups[0].QueryOffset, Is.EqualTo(0));
            Assert.That(result.PrecursorGroups[0].QueryCount, Is.EqualTo(2));
            Assert.That(result.PrecursorGroups[1].QueryOffset, Is.EqualTo(2));
            Assert.That(result.PrecursorGroups[1].QueryCount, Is.EqualTo(3));
        }

        [Test]
        public void Generate_PrecursorOutsideAllWindows_IsSkipped()
        {
            var p = new DiaSearchParameters { PpmTolerance = 20f, RtToleranceMinutes = 2f, MinFragmentsRequired = 1 };
            var precursors = new List<LibraryPrecursorInput>
            {
                MakePrecursor("OUTSIDE", 800.0, 2, new float[] { 200f }, new float[] { 100f }, rt: 11.0)
            };
            var result = DiaLibraryQueryGenerator.Generate(precursors, _scanIndex, p);
            Assert.That(result.Queries.Length, Is.EqualTo(0));
            Assert.That(result.SkippedNoWindow, Is.EqualTo(1));
        }

        [Test]
        public void Generate_NullRetentionTime_UsesFullRunRange()
        {
            var p = new DiaSearchParameters { PpmTolerance = 20f, RtToleranceMinutes = 2f, MinFragmentsRequired = 1 };
            var precursors = new List<LibraryPrecursorInput>
            {
                MakePrecursor("NORT", 412.5, 2, new float[] { 200f }, new float[] { 100f }, rt: null)
            };
            var result = DiaLibraryQueryGenerator.Generate(precursors, _scanIndex, p);
            Assert.That(result.Queries.Length, Is.EqualTo(1));
            Assert.That(result.Queries[0].RtMin, Is.EqualTo(_scanIndex.GetGlobalRtMin()).Within(0.01f));
            Assert.That(result.Queries[0].RtMax, Is.EqualTo(_scanIndex.GetGlobalRtMax()).Within(0.01f));
        }

        [Test]
        public void Generate_ZeroFragments_IsSkipped()
        {
            var p = new DiaSearchParameters { PpmTolerance = 20f, RtToleranceMinutes = 2f, MinFragmentsRequired = 1 };
            var precursors = new List<LibraryPrecursorInput>
            {
                MakePrecursor("NOFRAG", 412.5, 2, new float[0], new float[0], rt: 11.0)
            };
            var result = DiaLibraryQueryGenerator.Generate(precursors, _scanIndex, p);
            Assert.That(result.Queries.Length, Is.EqualTo(0));
            Assert.That(result.SkippedNoFragments, Is.EqualTo(1));
        }

        [Test]
        public void Generate_MixedValidAndSkipped()
        {
            var p = new DiaSearchParameters { PpmTolerance = 20f, RtToleranceMinutes = 2f, MinFragmentsRequired = 1 };
            var precursors = new List<LibraryPrecursorInput>
            {
                MakePrecursor("VALID", 412.5, 2, new float[] { 200f, 300f }, new float[] { 100f, 80f }, rt: 11.0),
                MakePrecursor("OUTSIDE", 800.0, 2, new float[] { 200f }, new float[] { 100f }, rt: 11.0),
                MakePrecursor("ALSOVALID", 475.0, 3, new float[] { 200f, 300f }, new float[] { 100f, 80f }, rt: 12.0)
            };
            var result = DiaLibraryQueryGenerator.Generate(precursors, _scanIndex, p);
            Assert.That(result.Queries.Length, Is.EqualTo(4));
            Assert.That(result.PrecursorGroups.Length, Is.EqualTo(2));
            Assert.That(result.SkippedNoWindow, Is.EqualTo(1));
            Assert.That(result.PrecursorGroups[0].InputIndex, Is.EqualTo(0));
            Assert.That(result.PrecursorGroups[1].InputIndex, Is.EqualTo(2));
        }

        [Test]
        public void Generate_EmptyList_ReturnsEmpty()
        {
            var p = new DiaSearchParameters { PpmTolerance = 20f, RtToleranceMinutes = 2f, MinFragmentsRequired = 1 };
            var result = DiaLibraryQueryGenerator.Generate(new List<LibraryPrecursorInput>(), _scanIndex, p);
            Assert.That(result.Queries.Length, Is.EqualTo(0));
            Assert.That(result.PrecursorGroups.Length, Is.EqualTo(0));
        }

        [Test]
        public void Generate_QueryIdsAreSequential()
        {
            var p = new DiaSearchParameters { PpmTolerance = 20f, RtToleranceMinutes = 2f, MinFragmentsRequired = 1 };
            var precursors = new List<LibraryPrecursorInput>
            {
                MakePrecursor("A", 412.5, 2, new float[] { 200f, 300f }, new float[] { 100f, 80f }, rt: 11.0),
                MakePrecursor("B", 437.5, 2, new float[] { 250f }, new float[] { 90f }, rt: 12.0)
            };
            var result = DiaLibraryQueryGenerator.Generate(precursors, _scanIndex, p);
            for (int i = 0; i < result.Queries.Length; i++)
                Assert.That(result.Queries[i].QueryId, Is.EqualTo(i));
        }

        // ═══════════════════════════════════════════════════════════════════
        // AssembleResults
        // ═══════════════════════════════════════════════════════════════════

        [Test]
        public void AssembleResults_PopulatesExtractedData()
        {
            var p = new DiaSearchParameters { PpmTolerance = 20f, RtToleranceMinutes = 2f, MinFragmentsRequired = 1 };
            var precursors = new List<LibraryPrecursorInput>
            {
                MakePrecursor("PEPTIDE", 412.5, 2, new float[] { 200f, 300f, 400f }, new float[] { 100f, 80f, 60f }, rt: 11.0)
            };
            var genResult = DiaLibraryQueryGenerator.Generate(precursors, _scanIndex, p);
            var extractionResults = new FragmentResult[]
            {
                new FragmentResult(0, 5, 0, 5, 500f),
                new FragmentResult(1, 3, 5, 3, 200f),
                new FragmentResult(2, 0, 8, 0, 0f),
            };
            var extractionResult = new ExtractionResult(extractionResults, Array.Empty<float>(), Array.Empty<float>(), 0);
            var results = DiaLibraryQueryGenerator.AssembleResults(precursors, genResult, extractionResult, p);
            Assert.That(results.Count, Is.EqualTo(1));
            Assert.That(results[0].Sequence, Is.EqualTo("PEPTIDE"));
            Assert.That(results[0].FragmentsDetected, Is.EqualTo(2));
            Assert.That(results[0].ExtractedIntensities[0], Is.EqualTo(500f));
            Assert.That(results[0].ExtractedIntensities[2], Is.EqualTo(0f));
        }

        [Test]
        public void AssembleResults_FiltersByMinFragments()
        {
            var p = new DiaSearchParameters { PpmTolerance = 20f, RtToleranceMinutes = 2f, MinFragmentsRequired = 3 };
            var precursors = new List<LibraryPrecursorInput>
            {
                MakePrecursor("PEPTIDE", 412.5, 2, new float[] { 200f, 300f, 400f }, new float[] { 100f, 80f, 60f }, rt: 11.0)
            };
            var genResult = DiaLibraryQueryGenerator.Generate(precursors, _scanIndex, p);
            var extractionResults = new FragmentResult[]
            {
                new FragmentResult(0, 5, 0, 5, 500f),
                new FragmentResult(1, 3, 5, 3, 200f),
                new FragmentResult(2, 0, 8, 0, 0f),
            };
            var extractionResult = new ExtractionResult(extractionResults, Array.Empty<float>(), Array.Empty<float>(), 0);
            var results = DiaLibraryQueryGenerator.AssembleResults(precursors, genResult, extractionResult, p);
            Assert.That(results.Count, Is.EqualTo(0));
        }

        [Test]
        public void AssembleResults_DecoyFlagPreserved()
        {
            var p = new DiaSearchParameters { PpmTolerance = 20f, RtToleranceMinutes = 2f, MinFragmentsRequired = 1 };
            var precursors = new List<LibraryPrecursorInput>
            {
                MakePrecursor("TARGET", 412.5, 2, new float[] { 200f }, new float[] { 100f }, rt: 11.0, isDecoy: false),
                MakePrecursor("DECOY", 437.5, 2, new float[] { 200f }, new float[] { 100f }, rt: 12.0, isDecoy: true)
            };
            var genResult = DiaLibraryQueryGenerator.Generate(precursors, _scanIndex, p);
            var extractionResults = new FragmentResult[]
            {
                new FragmentResult(0, 3, 0, 3, 100f),
                new FragmentResult(1, 3, 3, 3, 100f),
            };
            var extractionResult = new ExtractionResult(extractionResults, Array.Empty<float>(), Array.Empty<float>(), 0);
            var results = DiaLibraryQueryGenerator.AssembleResults(precursors, genResult, extractionResult, p);
            Assert.That(results.Count, Is.EqualTo(2));
            Assert.That(results[0].IsDecoy, Is.False);
            Assert.That(results[1].IsDecoy, Is.True);
        }

        // ═══════════════════════════════════════════════════════════════════
        // DiaSearchResult
        // ═══════════════════════════════════════════════════════════════════

        [Test]
        public void DiaSearchResult_Constructor_InitializesCorrectly()
        {
            var r = new DiaSearchResult("PEPTIDE", 2, 412.5, 0, false, 5, 12.0, 10f, 14f);
            Assert.That(r.Sequence, Is.EqualTo("PEPTIDE"));
            Assert.That(r.ChargeState, Is.EqualTo(2));
            Assert.That(r.ExtractedIntensities.Length, Is.EqualTo(5));
            Assert.That(float.IsNaN(r.DotProductScore), Is.True);
        }

        [Test]
        public void DiaSearchResult_MeetsMinFragments()
        {
            var r = new DiaSearchResult("PEP", 2, 400, 0, false, 5, 12.0, 10f, 14f);
            r.FragmentsDetected = 3;
            Assert.That(r.MeetsMinFragments(3), Is.True);
            Assert.That(r.MeetsMinFragments(4), Is.False);
        }

        [Test]
        public void DiaSearchResult_FragmentDetectionRate()
        {
            var r = new DiaSearchResult("PEP", 2, 400, 0, false, 4, 12.0, 10f, 14f);
            r.FragmentsDetected = 3;
            Assert.That(r.FragmentDetectionRate, Is.EqualTo(0.75f).Within(0.001f));
        }

        [Test]
        public void DiaSearchResult_ToString_IncludesKeyInfo()
        {
            var r = new DiaSearchResult("PEPTIDE", 2, 412.5, 0, true, 5, 12.0, 10f, 14f);
            r.DotProductScore = 0.85f;
            r.FragmentsDetected = 4;
            string s = r.ToString();
            Assert.That(s, Does.Contain("PEPTIDE/2"));
            Assert.That(s, Does.Contain("DECOY"));
            Assert.That(s, Does.Contain("4/5"));
        }

        // ═══════════════════════════════════════════════════════════════════
        // DiaSearchParameters
        // ═══════════════════════════════════════════════════════════════════

        [Test]
        public void DiaSearchParameters_Defaults()
        {
            var p = new DiaSearchParameters();
            Assert.That(p.PpmTolerance, Is.EqualTo(20f));
            Assert.That(p.RtToleranceMinutes, Is.EqualTo(5.0f));
            Assert.That(p.MinFragmentsRequired, Is.EqualTo(3));
            Assert.That(p.PreferGpu, Is.False);
        }

        [Test]
        public void DiaSearchParameters_EffectiveMaxThreads()
        {
            var p = new DiaSearchParameters { MaxThreads = -1 };
            Assert.That(p.EffectiveMaxThreads, Is.EqualTo(Environment.ProcessorCount));
            p.MaxThreads = 4;
            Assert.That(p.EffectiveMaxThreads, Is.EqualTo(4));
        }

        // ═══════════════════════════════════════════════════════════════════
        // LibraryPrecursorInput
        // ═══════════════════════════════════════════════════════════════════

        [Test]
        public void LibraryPrecursorInput_MismatchedArrays_Throws()
        {
            Assert.Throws<ArgumentException>(() =>
                new LibraryPrecursorInput("PEP", 400, 2, 12.0, false, new float[] { 200f, 300f }, new float[] { 100f }));
        }

        [Test]
        public void LibraryPrecursorInput_NullSequence_Throws()
        {
            Assert.Throws<ArgumentNullException>(() =>
                new LibraryPrecursorInput(null, 400, 2, 12.0, false, new float[] { 200f }, new float[] { 100f }));
        }

        [Test]
        public void LibraryPrecursorInput_FragmentCount()
        {
            var input = new LibraryPrecursorInput("PEP", 400, 2, 12.0, false,
                new float[] { 200f, 300f, 400f }, new float[] { 100f, 80f, 60f });
            Assert.That(input.FragmentCount, Is.EqualTo(3));
        }
    }
}