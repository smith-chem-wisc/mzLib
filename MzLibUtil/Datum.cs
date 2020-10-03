using System.Text;

namespace MzLibUtil
{
    public class Datum
    {
        public readonly double X;
        public readonly double? Y;
        public readonly double? Z;
        public readonly double? XError;
        public readonly double? YError;
        public readonly double? ZError;
        public readonly double Weight;
        public readonly string Label;

        public Datum(double x, double? y = null, double? z = null, double? xError = null, double? yError = null,
            double? zError = null, string label = null, double weight = 1)
        {
            this.X = x;
            this.Y = y;
            this.Z = z;

            this.XError = xError;
            this.YError = yError;
            this.ZError = zError;

            this.Label = label;
            this.Weight = weight;
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            sb.Append(X.ToString("F2"));

            if (XError.HasValue)
            {
                sb.Append(" ±");
                sb.Append(XError.Value);
            }

            if (Y.HasValue)
            {
                sb.Append(", ");
                sb.Append(Y.Value.ToString("F2"));
            }

            if (YError.HasValue)
            {
                sb.Append(" ±");
                sb.Append(YError.Value);
            }

            if (Z.HasValue)
            {
                sb.Append(", ");
                sb.Append(Z.Value.ToString("F2"));
            }

            if (ZError.HasValue)
            {
                sb.Append(" ±");
                sb.Append(ZError.Value);
            }

            if (Label != null)
            {
                sb.Append("; ");
                sb.Append(Label);
            }

            if (Weight != 1)
            {
                sb.Append("; ");
                sb.Append(Weight.ToString("F2"));
            }

            return sb.ToString();
        }
    }
}
