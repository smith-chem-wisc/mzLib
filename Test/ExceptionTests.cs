using NUnit.Framework;
using UsefulProteomicsDatabases;
using System;

namespace Test
{
    [TestFixture]
    public sealed class ExceptionTests
    {
        #region Public Methods

        [Test]
        public void PtmListLoaderExceptionTest1()
        {
            PtmListLoaderException p = new PtmListLoaderException();
            Assert.IsTrue(p is Exception);
        }

        [Test]
        public void PtmListLoaderExceptionTest2()
        {
            string message = "hey";
            Exception inner = new Exception();
            PtmListLoaderException p = new PtmListLoaderException(message, inner);
            Assert.IsTrue(p is Exception);
        }

        #endregion Public Methods
    }
}