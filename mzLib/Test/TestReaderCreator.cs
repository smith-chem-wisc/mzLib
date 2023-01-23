using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using MzLibUtil;
using NUnit.Framework;
using Readers; 
namespace Test;

[ExcludeFromCodeCoverage]
public class TestReaderCreator
{
    [Test]
    public void TestCreateReader()
    {
        // assert that all read successfully, and that an unknown file path returns null.
        string[] filePaths = new string[]
        {
            "small.raw", 
            "tester.mgf", 
            "SmallCalibratibleYeast.mzml", 
            "humanInsulin.fasta", 
            "fakeFile.d"
        };
        
        foreach (string path in filePaths)
        {
            string testPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", path);
            if (path == "fakeFile.d")
            {
                Assert.Throws<NotImplementedException>(() =>
                {
                    var reader = MsDataFileReader.CreateReader(testPath);
                }); 
                continue;
            }else if (path == "humanInsulin.fasta")
            {
                Assert.Throws<MzLibException>(() =>
                {
                    var reader = MsDataFileReader.CreateReader(testPath);
                });
                continue; 
            }
            var reader = MsDataFileReader.CreateReader(testPath);
            reader.LoadAllStaticData();
        }
        
    }
}