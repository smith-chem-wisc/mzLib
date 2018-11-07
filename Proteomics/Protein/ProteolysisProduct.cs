﻿namespace Proteomics
{
    public class ProteolysisProduct
    {
        public ProteolysisProduct(int? oneBasedBeginPosition, int? oneBasedEndPosition, string type)
        {
            OneBasedBeginPosition = oneBasedBeginPosition;
            OneBasedEndPosition = oneBasedEndPosition;
            Type = type ?? "";
        }

        public int? OneBasedBeginPosition { get; }
        public int? OneBasedEndPosition { get; }
        public string Type { get; }

        public override bool Equals(object obj)
        {
            ProteolysisProduct pp = obj as ProteolysisProduct;
            return pp != null
                && pp.OneBasedBeginPosition.Equals(OneBasedBeginPosition)
                && pp.OneBasedEndPosition.Equals(OneBasedEndPosition)
                && (pp.Type == null && Type == null || pp.Type.Equals(Type));
        }

        public override int GetHashCode()
        {
            return (OneBasedBeginPosition ?? 0).GetHashCode() 
                ^ (OneBasedEndPosition ?? 0).GetHashCode() 
                ^ Type.GetHashCode(); // null handled in constructor
        }
    }
}