using Easy.Common.Extensions;
using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Test.FileReadingTests
{
    public class TestMaxQuantResultFile
    {
        [TestFixture]
        [ExcludeFromCodeCoverage]
        public class MaxQuantEvidenceTests
        {
            [Test]
            public void TestMaxQuantEvidenceResultFile()
            {
                var path = @"D:\Kelly_TwoProteomeData\combined_MaxQuant_Kelly\txt\evidence.txt";
                var file = new MaxQuantEvidenceFile(path);
                file.LoadResults();

                var results = file.Results.ToList();

                Assert.That(results.IsNotNullOrEmpty(), "No results were loaded from the file");

            }
        }   
    }
}
