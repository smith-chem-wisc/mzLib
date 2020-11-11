using System.Text;

namespace MzLibUtil
{
    public class Datum
    {
        public readonly double X;
        public readonly double? Y;
        public readonly double? Z;
        public readonly double Weight;
        public readonly string Label;

        public Datum(double x, double? y = null, double? z = null, string label = null, double weight = 1)
        {
            this.X = x;
            this.Y = y;
            this.Z = z;

            this.Label = label;
            this.Weight = weight;
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            sb.Append(X);

            if (Y.HasValue)
            {
                sb.Append(", ");
                sb.Append(Y);
            }

            if (Z.HasValue)
            {
                sb.Append(", ");
                sb.Append(Z);
            }

            if (Label != null)
            {
                sb.Append("; ");
                sb.Append(Label);
            }

            if (Weight != 1)
            {
                sb.Append("; ");
                sb.Append(Weight);
            }

            return sb.ToString();
        }
    }
}
