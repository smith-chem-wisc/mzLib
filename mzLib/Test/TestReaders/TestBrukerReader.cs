using System;
using System.Diagnostics.CodeAnalysis;
using MathNet.Numerics.Distributions;
using NUnit.Framework;
using Readers; 
namespace Test;
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