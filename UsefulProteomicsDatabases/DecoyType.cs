using Proteomics.ProteolyticDigestion;
using System;

namespace UsefulProteomicsDatabases
{
    public enum DecoyType
    {
        None,
        Reverse,
        Slide,
        Shuffle,
        Random
    }

    public class DecoyTypeClass
    {
        public DecoyTypeClass(DecoyType decoyType, Protease protease)
        {
            DecoyType = decoyType;
            Protease = protease;
            if (DecoyType == DecoyType.Shuffle && protease == null)
            {
                throw new ArgumentException("DecoyType of shuffled must have a protease specified");
            }
        }

        public DecoyType DecoyType { get; private set; }
        public Protease Protease { get; private set; }
    }
}