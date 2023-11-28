using System;
using System.IO;
using System.Linq;
using MachineLearning;
using MachineLearning.RetentionTimePredictionModels;
using NUnit.Framework;
using UsefulProteomicsDatabases;

namespace Test
{
    public class TestChronologer
    {
        [Test]
        public void TestChronologerMethods()
        {
            var dictionary =
                ChronologerDictionary.GetChronologerDictionary(ChronologerDictionary.TypeOfDictionary.Chronologer);
            var chronologer = new Chronologer(dictionary, Path.Combine(AppDomain.CurrentDomain.BaseDirectory,
                "RetentionTimePredictionModels",
                "Chronologer_20220601193755_TorchSharp.dat"));
            
            //Default is Eval Mode
            Assert.That(chronologer.training == false);

            //Test TrainingMode 
            chronologer.TrainingMode();
            Assert.That(chronologer.training == true);
            
            //Test EvaluationMode
            chronologer.EvaluationMode();
            Assert.That(chronologer.training == false);

            //Check how many modules are in the model
            Assert.That(chronologer.modules().Count() == 24);

        }
    }
}
