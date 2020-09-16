namespace BayesianEstimation
{
    public class Datum
    {
        public readonly double Dimension1;
        public readonly double? Dimension2;
        public readonly double? Dimension3;
        public readonly double Weight;

        public Datum(double dataValue, double? dimension2 = null, double? dimension3 = null, double weight = 1)
        {
            this.Dimension1 = dataValue;
            this.Weight = weight;

            this.Dimension2 = dimension2;
            this.Dimension3 = dimension3;
        }
    }
}
