using NUnit.Framework;
using Proteomics;

namespace Test
{
    [TestFixture]
    public sealed class MyTests
    {
        [Test]
        public void ModificationCollectionTest()
        {
            OldSchoolModification mod1 = new OldSchoolModification(10, "mass 10 modification");
            OldSchoolModification mod2 = new OldSchoolModification(100, "mass 100 modification");
            OldSchoolModification mod3 = new OldSchoolModification(1000, "mass 1000 modification");
            ModificationCollection a = new ModificationCollection(mod1, mod2, mod3, mod1);
            ModificationCollection b = new ModificationCollection(mod1, mod3, mod1, mod2);
            Assert.IsTrue(a.Equals(b));
            ModificationCollection c = new ModificationCollection(mod1);
            Assert.IsFalse(c.Equals(b));
        }
    }
}