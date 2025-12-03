using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Test.Omics
{
    [ExcludeFromCodeCoverage]
    internal class IBioPolymerGroupTests
    {
        [Test]
        public void TestBioPolymerGroupEquality()
        {
            // Since IBioPolymerGroup is an interface, we need to create a mock or a simple implementation for testing
            var bioPolymerGroup1 = new MockBioPolymerGroup("GroupA");
            var bioPolymerGroup2 = new MockBioPolymerGroup("GroupA");
            var bioPolymerGroup3 = new MockBioPolymerGroup("GroupB");
            Assert.That(bioPolymerGroup1.Equals(bioPolymerGroup2), Is.True);
            Assert.That(bioPolymerGroup1.Equals(bioPolymerGroup3), Is.False);
            Assert.That(bioPolymerGroup1.Equals(null), Is.False);
        }
    }
}
