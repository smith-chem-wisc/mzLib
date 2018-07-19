using Chemistry;

namespace Proteomics.Fragmentation
{
    public class TheoreticalFragmentIon
    {
        /// <summary>
        /// Constructs a new TheoreticalFragmentIon given information about its theoretical properties
        /// </summary>
        public TheoreticalFragmentIon(double mass, double theorIntensity, int charge, ProductType productType, int ionNumber)
        {
            Mass = mass;
            Charge = charge;
            Intensity = theorIntensity;
            IonNumber = ionNumber;
            ProductType = productType;
            Mz = mass.ToMz(charge);
        }

        public double Mass { get; }
        public int Charge { get; }
        public double Mz { get; }
        public double Intensity { get; }
        public ProductType ProductType { get; }
        public int IonNumber { get; }

        /// <summary>
        /// Summarizes a TheoreticalFragmentIon into a string for debug purposes
        /// </summary>
        public override string ToString()
        {
            return ProductType + "" + IonNumber + ";" + Mass;
        }
    }
}
