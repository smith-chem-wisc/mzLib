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
            Modification mod1 = new Modification(10, "mass 10 modification");
            Modification mod2 = new Modification(100, "mass 100 modification");
            Modification mod3 = new Modification(1000, "mass 1000 modification");
            ModificationCollection a = new ModificationCollection(mod1, mod2, mod3, mod1);
            ModificationCollection b = new ModificationCollection(mod1, mod3, mod1, mod2);
            Assert.IsTrue(a.Equals(b));
            ModificationCollection c = new ModificationCollection(mod1);
            Assert.IsFalse(c.Equals(b));
        }
    }
}