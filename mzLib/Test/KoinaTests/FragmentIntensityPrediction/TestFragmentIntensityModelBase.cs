using NUnit.Framework;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.Util;
using System;
using System.Collections.Generic;
using System.ComponentModel;

namespace Test.KoinaTests
{
    /// <summary>
    /// Test harness for FragmentIntensityModel validation logic.
    /// Allows testing of base class validation without requiring a full model implementation.
    /// </summary>
    internal class TestFragmentIntensityModel : FragmentIntensityModel
    {
        public override string ModelName => "TestModel";
        public override int MaxBatchSize => 32;
        public override int MaxNumberOfBatchesPerRequest { get; init; } = 10;
        public override int ThrottlingDelayInMilliseconds { get; init; } = 100;
        public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 1000;
        public override int MaxPeptideLength => 30;
        public override int MinPeptideLength => 7;
        public override IncompatibleModHandlingMode ModHandlingMode { get; init; }
        public override IncompatibleParameterHandlingMode ParameterHandlingMode { get; init; }

        public override HashSet<int> AllowedPrecursorCharges { get; }
        public override HashSet<int> AllowedCollisionEnergies { get; }
        public override HashSet<string> AllowedInstrumentTypes { get; }
        public override HashSet<string> AllowedFragmentationTypes { get; }

        public TestFragmentIntensityModel(
            HashSet<int>? allowedCharges = null,
            HashSet<int>? allowedEnergies = null,
            HashSet<string>? allowedInstruments = null,
            HashSet<string>? allowedFragmentations = null,
            IncompatibleModHandlingMode modHandlingMode = IncompatibleModHandlingMode.ReturnNull,
            IncompatibleParameterHandlingMode parameterHandlingMode = IncompatibleParameterHandlingMode.ReturnNull
            )
        {
            AllowedPrecursorCharges = allowedCharges ?? new HashSet<int> { 2, 3, 4 };
            AllowedCollisionEnergies = allowedEnergies ?? new HashSet<int>();
            AllowedInstrumentTypes = allowedInstruments ?? new HashSet<string>();
            AllowedFragmentationTypes = allowedFragmentations ?? new HashSet<string>();
            ModHandlingMode = modHandlingMode;
            ParameterHandlingMode = parameterHandlingMode;
        }

        protected override List<Dictionary<string, object>> ToBatchedRequests(List<FragmentIntensityPredictionInput> validInputs)
        {
            // Minimal implementation for testing
            return new List<Dictionary<string, object>>();
        }

        // Expose protected method for direct testing
        public bool TestValidateModelSpecificInputs(FragmentIntensityPredictionInput input, out WarningException? warning)
        {
            return ValidateModelSpecificInputs(input, out warning);
        }
    }

    public class FragmentIntensityModelValidationTests
    {
        [Test]
        public void TestCollisionEnergyValidation()
        {
            var allowedEnergies = new HashSet<int> { 25, 30, 35 };

            var validInput = new FragmentIntensityPredictionInput("PEPTIDEK", 2, 30, null, null);
            var invalidInput = new FragmentIntensityPredictionInput("PEPTIDEK", 2, 40, null, null);
            var nullInput = new FragmentIntensityPredictionInput("PEPTIDEK", 2, null, null, null);

            var model = new TestFragmentIntensityModel(allowedEnergies: allowedEnergies);
            var isValid = model.TestValidateModelSpecificInputs(validInput, out var validWarning);
            var isInvalid = model.TestValidateModelSpecificInputs(invalidInput, out var invalidWarning);
            var isNull = model.TestValidateModelSpecificInputs(nullInput, out var nullWarning); // null energy should throw because constraints are defined.

            Assert.That(isValid, Is.True);
            Assert.That(isInvalid, Is.False);
            Assert.That(validWarning, Is.Null);
            Assert.That(invalidWarning, Is.Not.Null);
            Assert.That(isNull, Is.False);
            Assert.That(nullWarning, Is.Not.Null);

            // Test behavior when handling mode is set to throw exception
            model = new TestFragmentIntensityModel(allowedEnergies: allowedEnergies, parameterHandlingMode: IncompatibleParameterHandlingMode.ThrowException);
            Assert.DoesNotThrow(() => model.TestValidateModelSpecificInputs(validInput, out _));
            Assert.Throws<ArgumentException>(() => model.TestValidateModelSpecificInputs(invalidInput, out _));
            Assert.Throws<ArgumentException>(() => model.TestValidateModelSpecificInputs(nullInput, out _));

            // Test null behavior when CE constraint not specified
            model = new TestFragmentIntensityModel(parameterHandlingMode: IncompatibleParameterHandlingMode.ReturnNull);
            Assert.DoesNotThrow(() => model.TestValidateModelSpecificInputs(nullInput, out var warning));
            model = new TestFragmentIntensityModel(parameterHandlingMode: IncompatibleParameterHandlingMode.ThrowException);
            Assert.DoesNotThrow(() => model.TestValidateModelSpecificInputs(nullInput, out var warningWithConstraints));
        }

        [Test]
        public void TestInstrumentTypeValidation()
        {
            var allowedInstruments = new HashSet<string> { "QE", "LUMOS", "VELOS" };

            var validInput = new FragmentIntensityPredictionInput("PEPTIDEK", 2, null, "QE", null);
            var invalidInput = new FragmentIntensityPredictionInput("PEPTIDEK", 2, null, "UNKNOWN", null);
            var nullInput = new FragmentIntensityPredictionInput("PEPTIDEK", 2, null, null, null);

            var model = new TestFragmentIntensityModel(allowedInstruments: allowedInstruments);
            var isValid = model.TestValidateModelSpecificInputs(validInput, out var validWarning);
            var isInvalid = model.TestValidateModelSpecificInputs(invalidInput, out var invalidWarning);
            var isNull = model.TestValidateModelSpecificInputs(nullInput, out var nullWarning); // null instrument should fail because constraints are defined.

            Assert.That(isValid, Is.True);
            Assert.That(isInvalid, Is.False);
            Assert.That(validWarning, Is.Null);
            Assert.That(invalidWarning, Is.Not.Null);
            Assert.That(invalidWarning.Message, Does.Contain("Instrument type 'UNKNOWN' is not supported"));
            Assert.That(isNull, Is.False);
            Assert.That(nullWarning, Is.Not.Null);

            // Test behavior when handling mode is set to throw exception
            model = new TestFragmentIntensityModel(allowedInstruments: allowedInstruments, parameterHandlingMode: IncompatibleParameterHandlingMode.ThrowException);
            Assert.DoesNotThrow(() => model.TestValidateModelSpecificInputs(validInput, out _));
            Assert.Throws<ArgumentException>(() => model.TestValidateModelSpecificInputs(invalidInput, out _));
            Assert.Throws<ArgumentException>(() => model.TestValidateModelSpecificInputs(nullInput, out _));

            // Test null behavior when instrument constraint not specified
            model = new TestFragmentIntensityModel(parameterHandlingMode: IncompatibleParameterHandlingMode.ReturnNull);
            Assert.DoesNotThrow(() => model.TestValidateModelSpecificInputs(nullInput, out var warning));
            model = new TestFragmentIntensityModel(parameterHandlingMode: IncompatibleParameterHandlingMode.ThrowException);
            Assert.DoesNotThrow(() => model.TestValidateModelSpecificInputs(nullInput, out var warningWithConstraints));
        }

        [Test]
        public void TestFragmentationTypeValidation()
        {
            var allowedFragmentations = new HashSet<string> { "HCD", "CID", "ETD" };

            var validInput = new FragmentIntensityPredictionInput("PEPTIDEK", 2, null, null, "HCD");
            var invalidInput = new FragmentIntensityPredictionInput("PEPTIDEK", 2, null, null, "ECD");
            var nullInput = new FragmentIntensityPredictionInput("PEPTIDEK", 2, null, null, null);

            var model = new TestFragmentIntensityModel(allowedFragmentations: allowedFragmentations);
            var isValid = model.TestValidateModelSpecificInputs(validInput, out var validWarning);
            var isInvalid = model.TestValidateModelSpecificInputs(invalidInput, out var invalidWarning);
            var isNull = model.TestValidateModelSpecificInputs(nullInput, out var nullWarning); // null fragmentation should fail because constraints are defined.

            Assert.That(isValid, Is.True);
            Assert.That(isInvalid, Is.False);
            Assert.That(validWarning, Is.Null);
            Assert.That(invalidWarning, Is.Not.Null);
            Assert.That(invalidWarning.Message, Does.Contain("Fragmentation type 'ECD' is not supported"));
            Assert.That(isNull, Is.False);
            Assert.That(nullWarning, Is.Not.Null);

            // Test behavior when handling mode is set to throw exception
            model = new TestFragmentIntensityModel(allowedFragmentations: allowedFragmentations, parameterHandlingMode: IncompatibleParameterHandlingMode.ThrowException);
            Assert.DoesNotThrow(() => model.TestValidateModelSpecificInputs(validInput, out _));
            Assert.Throws<ArgumentException>(() => model.TestValidateModelSpecificInputs(invalidInput, out _));
            Assert.Throws<ArgumentException>(() => model.TestValidateModelSpecificInputs(nullInput, out _));

            // Test null behavior when fragmentation constraint not specified
            model = new TestFragmentIntensityModel(parameterHandlingMode: IncompatibleParameterHandlingMode.ReturnNull);
            Assert.DoesNotThrow(() => model.TestValidateModelSpecificInputs(nullInput, out var warning));
            model = new TestFragmentIntensityModel(parameterHandlingMode: IncompatibleParameterHandlingMode.ThrowException);
            Assert.DoesNotThrow(() => model.TestValidateModelSpecificInputs(nullInput, out var warningWithConstraints));
        }

        [Test]
        public void TestPrecursorChargeValidation()
        {
            var allowedCharges = new HashSet<int> { 2, 3, 4 };

            var validInput = new FragmentIntensityPredictionInput("PEPTIDEK", 3, null, null, null);
            var invalidInput = new FragmentIntensityPredictionInput("PEPTIDEK", 5, null, null, null);

            var model = new TestFragmentIntensityModel(allowedCharges: allowedCharges);
            var isValid = model.TestValidateModelSpecificInputs(validInput, out var validWarning);
            var isInvalid = model.TestValidateModelSpecificInputs(invalidInput, out var invalidWarning);

            Assert.That(isValid, Is.True);
            Assert.That(isInvalid, Is.False);
            Assert.That(validWarning, Is.Null);
            Assert.That(invalidWarning, Is.Not.Null);
            Assert.That(invalidWarning.Message, Does.Contain("Precursor charge"));

            // Test behavior when handling mode is set to throw exception
            model = new TestFragmentIntensityModel(allowedCharges: allowedCharges, parameterHandlingMode: IncompatibleParameterHandlingMode.ThrowException);
            Assert.DoesNotThrow(() => model.TestValidateModelSpecificInputs(validInput, out _));
            Assert.Throws<ArgumentException>(() => model.TestValidateModelSpecificInputs(invalidInput, out _));
        }

        [Test]
        public void TestMultipleConstraints_AllValid_ReturnsTrue()
        {
            var model = new TestFragmentIntensityModel(
                allowedCharges: new HashSet<int> { 2, 3 },
                allowedEnergies: new HashSet<int> { 30, 35 },
                allowedInstruments: new HashSet<string> { "QE" },
                allowedFragmentations: new HashSet<string> { "HCD" }
            );

            var input = new FragmentIntensityPredictionInput("PEPTIDEK", 2, 30, "QE", "HCD");

            var isValid = model.TestValidateModelSpecificInputs(input, out var warning);

            Assert.That(isValid, Is.True);
            Assert.That(warning, Is.Null);
        }

        [Test]
        public void TestMultipleConstraints_OneInvalid_ReturnsFalse()
        {
            var model = new TestFragmentIntensityModel(
                allowedCharges: new HashSet<int> { 2, 3 },
                allowedEnergies: new HashSet<int> { 30, 35 },
                allowedInstruments: new HashSet<string> { "QE" },
                allowedFragmentations: new HashSet<string> { "HCD" }
            );

            var input = new FragmentIntensityPredictionInput("PEPTIDEK", 2, 40, "QE", "HCD");

            var isValid = model.TestValidateModelSpecificInputs(input, out var warning);

            Assert.That(isValid, Is.False);
            Assert.That(warning, Is.Not.Null);
            Assert.That(warning.Message, Does.Contain("Collision energy"));
        }

        [Test]
        public void TestValidation_EmptyConstraintSets_AllowsAnyValue()
        {
            // When constraint sets are empty, validation should pass (no restrictions)
            var model = new TestFragmentIntensityModel(
                allowedEnergies: new HashSet<int>(),
                allowedInstruments: new HashSet<string>(),
                allowedFragmentations: new HashSet<string>()
            );

            var input = new FragmentIntensityPredictionInput("PEPTIDEK", 2, 999, "AnyInstrument", "AnyFragType");

            var isValid = model.TestValidateModelSpecificInputs(input, out var warning);

            Assert.That(isValid, Is.True);
            Assert.That(warning, Is.Null);
        }
    }
}