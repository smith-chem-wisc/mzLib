using NetSerializer;
using Omics.Digestion;
using Omics.Modifications;

namespace Omics;

/// <summary>
/// An interface to represent the required fields for serialization of a Full Sequence containing object in MetaMorpheus
/// </summary>
public interface ISerializableSequence
{
    Type[] GetTypesToSerialize();
    void SetNonSerializedPeptideInfo(IDictionary<string, Modification> allKnownMods, IDictionary<string, IBioPolymer> accessionToProtein, IDigestionParams digestionParams); 
}

public static class SerializableSequenceExtensions
{
    public static Serializer GetSequenceSerializer<T>(this T toSerialize) where T : ISerializableSequence
    {
        var collectionType = typeof(List<T>);
        var implementationSpecificTypes = toSerialize.GetTypesToSerialize();

        var types = new List<Type>(implementationSpecificTypes) { collectionType };
        var serializer = new Serializer(types);
        return serializer;
    }

    public static Serializer GetSequenceSerializer<T>(this List<T> toSerialize) where T : ISerializableSequence
    {
        T instance = toSerialize.Count == 0 ?
            (T)Activator.CreateInstance(typeof(T), true)!
            : toSerialize[0];

        var collectionType = toSerialize.GetType();
        var implementationSpecificTypes = instance.GetTypesToSerialize();

        var types = new List<Type>(implementationSpecificTypes) { collectionType };
        var serializer = new Serializer(types);
        return serializer;
    }
}
