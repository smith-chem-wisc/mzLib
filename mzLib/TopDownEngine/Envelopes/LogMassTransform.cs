using Chemistry;
using System;

namespace TopDownEngine.Envelopes;

public static class LogMassTransform
{
    public static double LogMz(double mz) => Math.Log(mz - Constants.ProtonMass);

    public static double[] BuildTemplate(int zMax)
    {
        double[] template = new double[zMax];
        for (int z = 1; z <= zMax; z++)
        {
            template[z - 1] = -Math.Log(z);
        }

        return template;
    }
}
