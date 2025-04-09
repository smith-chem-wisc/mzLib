using NetSerializer;
using Omics.Digestion;
using Omics.Modifications;

namespace Omics;

/// <summary>
/// An interface to represent the required fields for serialization of a Full Sequence containing object in MetaMorpheus
/// </summary>
public interface ISerializableSequence
{
    /// <summary>
    /// Types required for Nett Serializer to function properly
    /// </summary>
    /// <returns>An array of types required for serialization.</returns>
    Type[] GetTypesToSerialize();

    /// <summary>
    /// Set non-serialized values that are retained in MetaMorpheus during runtime
    /// </summary>
    /// <param name="allKnownMods">Dictionary of all known modifications.</param>
    /// <param name="accessionToProtein">Dictionary mapping accession to protein.</param>
    /// <param name="digestionParams">Parameters for digestion.</param>
    void SetNonSerializedPeptideInfo(IDictionary<string, Modification> allKnownMods, IDictionary<string, IBioPolymer> accessionToProtein, IDigestionParams digestionParams);

    /// <summary>
    /// Gets the types required for serialization for a specific implementation of ISerializableSequence.
    /// </summary>
    /// <typeparam name="T">The type implementing ISerializableSequence.</typeparam>
    /// <returns>An array of types required for serialization.</returns>
    static Type[] GetTypesToSerialize<T>() where T : ISerializableSequence
    {
        // Create an instance of the type T using reflection and a private empty constructor. 
        var instance = Activator.CreateInstance(typeof(T), true);
        if (instance == null)
        {
            throw new InvalidOperationException($"Cannot create an instance of {typeof(T)}. No Empty Constructor");
        }

        // Cast the instance to ISerializableSequence and call GetTypesToSerialize.
        var castedInstance = (ISerializableSequence)instance;
        return castedInstance.GetTypesToSerialize();
    }

    /// <summary>
    /// Creates a serializer for a specific implementation of ISerializableSequence.
    /// </summary>
    /// <typeparam name="T">The type implementing ISerializableSequence.</typeparam>
    /// <returns>A serializer for the specified type.</returns>
    static Serializer GetSequenceSerializer<T>() where T : ISerializableSequence
    {
        var types = GetTypesToSerialize<T>();
        var serializer = new Serializer(types);
        return serializer;
    }

    /// <summary>
    /// Creates a serializer for a list of a specific implementation of ISerializableSequence.
    /// </summary>
    /// <typeparam name="T">The type implementing ISerializableSequence.</typeparam>
    /// <param name="values">A list of values to be serialized.</param>
    /// <returns>A serializer for the specified type.</returns>
    static Serializer GetSequenceSerializer<T>(List<T> values) where T : ISerializableSequence
    {
        var types = GetTypesToSerialize<T>();
        var serializer = new Serializer(types);
        return serializer;
    }
}
