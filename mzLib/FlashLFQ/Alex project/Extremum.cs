using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ.Alex_project
{

    public enum ExtremumType { Minimum, Maximum }; //The spectra file the XIC came from

    internal class Extremum : IComparable<Extremum>, IEquatable<Extremum>
    {
        public readonly double Intensity;     //The aligned intensity of the Extremum point
        public readonly double RetentionTime; //The interpolated intensity of the Extremum point
        public readonly ExtremumType Type;    //The type of the Extremum point

        public Extremum(double intensity, double retentionTime, ExtremumType type)
        {
            Intensity = intensity;
            RetentionTime = retentionTime;
            Type = type;
        }



        
        public int CompareTo(Extremum others)
        {
            if (others == null)
            {
                return 1;
            }

            if (this - others > 0)
            {
                return 1;
            }

            else if (this - others < 0)
            {
                return -1;
            }

            else 
            {
                return 0;
            }

        }

        public bool Equals(Extremum others)
        {
            if (this.Intensity == others.Intensity && Math.Abs(this - others) < 0.006 && this.Type == others.Type)
            {
                return true;
            }

            return false;
        }

        public override bool Equals(Object obj)
        {
            return Equals(this, (Extremum)obj);
        }

        public static double operator - (Extremum extremun1, Extremum extrenum2)
        {
            double rtDiff = extremun1.RetentionTime - extrenum2.RetentionTime;
            return rtDiff;
        }

    
    }
}
