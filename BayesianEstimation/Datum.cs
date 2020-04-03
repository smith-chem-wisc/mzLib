namespace BayesianEstimation
{
    public class Datum
    {
        public readonly double DataValue;
        public readonly double Weight;

        public Datum(double dataValue, double weight)
        {
            this.DataValue = dataValue;
            this.Weight = weight;
        }
    }
}
