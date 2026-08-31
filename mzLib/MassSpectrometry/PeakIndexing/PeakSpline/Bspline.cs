using Easy.Common.Extensions;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static System.Runtime.InteropServices.JavaScript.JSType;
using MzLibUtil;

namespace MassSpectrometry
{
    /// <summary>
    /// This Bspline class is copied from DIA-Umpire code. It is a more generalized spline that allows different smoothing degrees.
    /// In theory, it could replace linear and cubic spline where the smoothing degree is explicitly set to 1 and 3 respectively.
    /// Only used for RT-based spline for now.
    /// </summary>
    public class Bspline : XicSpline
    {
        public int SmoothDegree { get; set; }
        public int NumberOfPoints { get; set; }

        public Bspline(int smoothDegree, int numberOfPoints)
        {
            SmoothDegree = smoothDegree;
            NumberOfPoints = numberOfPoints;
        }

        public override (double, double)[] GetXicSplineData(float[] rtArray, float[] intensityArray, double start = -1, double end = -1)
        {
            List<(double, double)> bsplineCollection = new List<(double, double)>();
            int p = SmoothDegree;
            int n = rtArray.Length - 1;
            int m = rtArray.Length + p;

            if (rtArray.Length <= p)
            {
                throw new MzLibException("The number of points in the input array must be greater than the degree of the Bspline.");
            }

            float[] bspline_T_ = new float[m + p];
            for (int i = 0; i <= n; i++)
            {
                bspline_T_[i] = 0;
                bspline_T_[m - i] = 1;
            }
            float intv = 1.0f / (m - 2 * p);
            for (int i = 1; i <= m - 1; i++)
            {
                bspline_T_[p + i] = bspline_T_[p + i - 1] + intv;
            }
            for (int i = 0; i < NumberOfPoints; i++)
            {
                float t = (float)i / NumberOfPoints;
                var pt = getbspline(rtArray, intensityArray, t, n, p, bspline_T_);
                bsplineCollection.Add(pt);
            }
            if (bsplineCollection[bsplineCollection.Count() - 1].Item1 < rtArray[rtArray.Length - 1])
            {
                bsplineCollection.Add((rtArray[rtArray.Length - 1], intensityArray[intensityArray.Length - 1]));
            }
            if (bsplineCollection[0].Item1 > rtArray[0])
            {
                bsplineCollection.Add((rtArray[0], intensityArray[0]));
            }
            return bsplineCollection.ToArray();
        }

        public (float, float) getbspline(float[] rtArray, float[] intensityArray, float t, int n, int p, float[] bspline_T_)
        {
            float x = 0, y = 0;
            for (int i = 0; i <= n; i++)
            {
                float a = bspline_base(i, p, t, bspline_T_);
                x += rtArray[i] * a;
                y += intensityArray[i] * a;
            }
            return new(x, y);
        }

        public float bspline_base(int i, int p, float t, float[] bspline_T_)
        {
            float n, c1, c2;
            float tn1 = 0;
            float tn2 = 0;
            if (p == 0)
            {
                if (bspline_T_[i] <= t && t < bspline_T_[i + 1] && bspline_T_[i] < bspline_T_[i + 1])
                {
                    n = 1;
                }
                else
                {
                    n = 0;
                }
            }
            else
            {
                if (bspline_T_[i + p] - bspline_T_[i] == 0)
                {
                    c1 = 0;
                }
                else
                {
                    tn1 = bspline_base(i, p - 1, t, bspline_T_);
                    c1 = (t - bspline_T_[i]) / (bspline_T_[i + p] - bspline_T_[i]);
                }
                if (bspline_T_[i + p + 1] - bspline_T_[i + 1] == 0)
                {
                    c2 = 0;
                }
                else
                {
                    tn2 = bspline_base(i + 1, p - 1, t, bspline_T_);
                    c2 = (bspline_T_[i + p + 1] - t) / (bspline_T_[i + p + 1] - bspline_T_[i + 1]);
                }
                n = c1 * tn1 + c2 * tn2;
            }
            return n;
        }
    }
}
