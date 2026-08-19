using System;
using System.Linq;
using System.Reflection;
using mzIdentML111.Generated;
using NUnit.Framework;

namespace Test.MzIdentML
{
    /// <summary>
    /// The mzIdentML reader ships three structurally-identical, XSD-generated DTO trees
    /// (namespaces mzIdentML110/111/120.Generated). Hand-writing per-property set/get tests
    /// for every version is what historically tanked coverage on these files. Instead, a single
    /// reflection-driven test walks every concrete generated type in a namespace and exercises
    /// each property's getter and setter. Concrete subclasses inherit the abstract base accessors
    /// (AbstractContactType, IdentifiableType, AbstractParamType, ...), so those are covered too.
    /// One parameterized test therefore covers all three versions.
    /// </summary>
    [TestFixture]
    public class MzIdentMlReflectionCoverageTests
    {
        [Test]
        [TestCase("mzIdentML110.Generated")]
        [TestCase("mzIdentML111.Generated")]
        [TestCase("mzIdentML120.Generated")]
        public void EveryGeneratedType_PropertyGetSet_RoundTrips(string generatedNamespace)
        {
            Assembly assembly = typeof(MzIdentMLType111).Assembly;

            Type[] concreteTypes = assembly.GetTypes()
                .Where(t => t.Namespace == generatedNamespace
                            && t.IsClass
                            && !t.IsAbstract
                            && t.GetConstructor(Type.EmptyTypes) != null)
                .ToArray();

            Assert.That(concreteTypes.Length, Is.GreaterThan(0),
                $"No instantiable generated types found in {generatedNamespace}.");

            foreach (Type type in concreteTypes)
            {
                object instance = Activator.CreateInstance(type);

                foreach (PropertyInfo property in type.GetProperties(BindingFlags.Public | BindingFlags.Instance))
                {
                    if (!property.CanRead || !property.CanWrite || property.GetIndexParameters().Length > 0)
                        continue;

                    object value = CreateValue(property.PropertyType);

                    // A non-nullable value type can never be assigned null; CreateValue always
                    // yields a non-null for those, so a null here means an unsupported reference type.
                    property.SetValue(instance, value);          // exercises the setter
                    object readBack = property.GetValue(instance); // exercises the getter

                    Assert.That(readBack, Is.EqualTo(value),
                        $"{type.Name}.{property.Name} did not round-trip.");
                }
            }
        }

        /// <summary>
        /// Produces a representative value for a property type so the setter/getter can be exercised.
        /// Abstract/interface reference types that can't be instantiated fall back to null, which still
        /// exercises the accessors.
        /// </summary>
        private static object CreateValue(Type type)
        {
            Type underlying = Nullable.GetUnderlyingType(type);
            if (underlying != null)
                return CreateValue(underlying);

            if (type == typeof(string))
                return "x";
            if (type == typeof(bool))
                return true;
            if (type == typeof(DateTime))
                return new DateTime(2020, 1, 1);
            if (type.IsEnum)
                return Enum.GetValues(type).GetValue(0);
            if (type == typeof(int) || type == typeof(long) || type == typeof(short) || type == typeof(byte) ||
                type == typeof(uint) || type == typeof(ulong) || type == typeof(ushort) || type == typeof(sbyte))
                return Convert.ChangeType(1, type);
            if (type == typeof(double) || type == typeof(float) || type == typeof(decimal))
                return Convert.ChangeType(1, type);
            if (type.IsValueType)
                return Activator.CreateInstance(type); // any other struct -> default
            if (type.IsArray)
                return Array.CreateInstance(type.GetElementType(), 0);
            if (type == typeof(object))
                return "x";

            // Reference type: instantiate when there is an accessible parameterless ctor.
            if (!type.IsAbstract && !type.IsInterface && type.GetConstructor(Type.EmptyTypes) != null)
            {
                try { return Activator.CreateInstance(type); }
                catch { return null; }
            }

            return null; // abstract/interface/no-ctor -> leave null; accessors still run
        }
    }
}
