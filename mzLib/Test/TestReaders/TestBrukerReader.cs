using System;
using System.Diagnostics.CodeAnalysis;
using NUnit.Framework;
using Readers;
using Readers.ReaderFactories;

namespace Test.TestReaders;
[ExcludeFromCodeCoverage]
public class TestBrukerReader
{
    [Test]
    [ExcludeFromCodeCoverage]
    public void TestBrukerFactory()
    {
        string dummyPath = "fakeFilePath.d";
        Assert.Throws<NotImplementedException>(() =>
        {
            BrukerReaderFactory bruk = new BrukerReaderFactory(dummyPath);
        }); 

    }
}