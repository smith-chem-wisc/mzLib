using Easy.Common.Extensions;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TorchSharp;
using UsefulProteomicsDatabases;
using Chronologer = MachineLearning.RetentionTimePredictionModels.Chronologer;

namespace Test
{
    public class TestChronologer
    {
        [Test]
        public void TestChronologerMethods()
        {
            var dictionary =
                DictionaryBuilder.GetChronologerDictionary(DictionaryBuilder.TypeOfDictionary.Chronologer);
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

            //Test forward Method
            Assert.That(chronologer.forward(torch.zeros(1, 52, torch.ScalarType.Int64)).GetType() == typeof(torch.Tensor));

            //Test ModificationDictionary getter
            Assert.That(chronologer.ModificationDictionary.IsNotNullOrEmpty());

        }

        [Test]
        public void TestDatasets()
        {
            var dictionary =
                DictionaryBuilder.GetChronologerDictionary(DictionaryBuilder.TypeOfDictionary.Chronologer);

            var chronologer = new Chronologer(dictionary, Path.Combine(AppDomain.CurrentDomain.BaseDirectory,
                "RetentionTimePredictionModels",
                "Chronologer_20220601193755_TorchSharp.dat"));

            List<string> warnings = new List<string>();

            var psms = Readers.SpectrumMatchTsvReader.ReadPsmTsv(
                Path.Combine(TestContext.CurrentContext.TestDirectory,
                    @"MachineLearningTests/AllPeptides.psmtsv"), out warnings);

            //Make sure the dataset is not null
            chronologer.CreateDataSet(psms, 0.1f, 0.1f, 64);

            //Test TrainDataset getter
            Assert.That(chronologer.TrainingDataset.Count > 0);

            //Test TrainDataset setter
            chronologer.TrainingDataset = null;
            Assert.That(chronologer.TrainingDataset == null);

            //Test TestingDataset getter
            Assert.That(chronologer.TestingDataset.Count > 0);

            //Test TestingDataset setter
            chronologer.TestingDataset = null;
            Assert.That(chronologer.TestingDataset == null);

            //Test ValidationDataset getter
            Assert.That(chronologer.ValidationDataset.Count > 0);

            //Test ValidationDataset setter
            chronologer.ValidationDataset = null;
            Assert.That(chronologer.ValidationDataset == null);

            //Test TrainDataLoader getter
            Assert.That(chronologer.TrainDataLoader.IsNotNullOrEmpty());

            //Test TrainDataLoader setter
            chronologer.TrainDataLoader = null;
            Assert.That(chronologer.TrainDataLoader == null);

            //Test TestDataLoader getter
            Assert.That(chronologer.TestDataLoader.IsNotNullOrEmpty());

            //Test TestDataLoader setter
            chronologer.TestDataLoader = null;
            Assert.That(chronologer.TestDataLoader == null);

            //Test ValidationDataLoader getter
            Assert.That(chronologer.ValidationDataLoader.IsNotNullOrEmpty());

            //Test ValidationDataLoader setter
            chronologer.ValidationDataLoader = null;
            Assert.That(chronologer.ValidationDataLoader == null);

            //Test Chronologer Constructor TrainingMode
            var chronologer2 = new Chronologer(dictionary, Path.Combine(AppDomain.CurrentDomain.BaseDirectory,
                               "RetentionTimePredictionModels",
                                              "Chronologer_20220601193755_TorchSharp.dat"), false);
            Assert.That(chronologer2.training == true);
        }

        [Test]
        public void TestTensorize()
        {
            var dictionary =
                DictionaryBuilder.GetChronologerDictionary(DictionaryBuilder.TypeOfDictionary.Chronologer);

            var chronologer = new Chronologer(dictionary, Path.Combine(AppDomain.CurrentDomain.BaseDirectory,
                               "RetentionTimePredictionModels",
                                              "Chronologer_20220601193755_TorchSharp.dat"));

            var psms = Readers.SpectrumMatchTsvReader.ReadPsmTsv(
                               Path.Combine(TestContext.CurrentContext.TestDirectory,
                                   @"MachineLearningTests/AllPeptides.psmtsv"), out List<string> warnings);

            chronologer.CreateDataSet(psms, 0.1f, 0.1f, 64);

            //Test Tensorize
            Assert.Throws<ArgumentException>(() => chronologer.Tensorize(1));

        }

        [Test]
        public void TestForwardPass()
        {
            var dictionary =
                DictionaryBuilder.GetChronologerDictionary(DictionaryBuilder.TypeOfDictionary.Chronologer);

            var chronologer = new Chronologer(dictionary, Path.Combine(AppDomain.CurrentDomain.BaseDirectory,
                "RetentionTimePredictionModels",
                "Chronologer_20220601193755_TorchSharp.dat"));

            var psms = Readers.SpectrumMatchTsvReader.ReadPsmTsv(
                Path.Combine(TestContext.CurrentContext.TestDirectory,
                    @"MachineLearningTests/AllPeptides.psmtsv"), out List<string> warnings);

            var prediction = chronologer.forward(torch.zeros(1, 52, torch.ScalarType.Int64));

            Assert.That(prediction.GetType() == typeof(torch.Tensor));
        }

        [Test]
        public void TestTrainMethod()
        {
            var dictionary =
                DictionaryBuilder.GetChronologerDictionary(DictionaryBuilder.TypeOfDictionary.Chronologer);

            var chronologer = new Chronologer(dictionary, Path.Combine(AppDomain.CurrentDomain.BaseDirectory,
                               "RetentionTimePredictionModels",
                                              "Chronologer_20220601193755_TorchSharp.dat"));

            var psms = Readers.SpectrumMatchTsvReader.ReadPsmTsv(
                               Path.Combine(TestContext.CurrentContext.TestDirectory,
                                   @"MachineLearningTests/AllPeptides.psmtsv"), out List<string> warnings);

            chronologer.CreateDataSet(psms, 0.1f, 0.1f, 64);

            //Test the case where the path does not contain .dat extension
            chronologer.Train("testingTrainingMethod", psms, dictionary, DeviceType.CPU, 0.1f,
                0.1f, 64, 50, 5);

            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, "testingTrainingMethod.dat"));
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, "testingTrainingMethod.txt"));

            //Test the case where the path has a .dat extension
            chronologer.Train("testingTrainingMethod.dat", psms, dictionary, DeviceType.CPU, 0.1f,
                0.1f, 64, 50, 5);

            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, "testingTrainingMethod.dat"));
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, "testingTrainingMethod.txt"));
        }
    }
}
