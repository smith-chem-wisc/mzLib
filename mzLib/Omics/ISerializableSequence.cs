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
        var types = toSerialize.GetTypesToSerialize();
        var serializer = new Serializer(types);
        return serializer;
    }
}
