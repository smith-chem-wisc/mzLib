

namespace MassSpectrometry.Dia.Calibration
{
    /// <summary>
    /// Wraps a LowessRtModel but overrides SigmaMinutes for window generation.
    /// Preserves the nonlinear ToMinutes() curve shape from the bootstrap.
    /// </summary>
    internal sealed class LowessSigmaOverrideModel : IRtCalibrationModel
    {
        private readonly LowessRtModel _inner;

        public LowessSigmaOverrideModel(LowessRtModel inner, double sigmaOverride)
        {
            _inner = inner;
            SigmaMinutes = sigmaOverride;
        }

        public double SigmaMinutes { get; }
        public double RSquared => _inner.RSquared;
        public double Slope => _inner.Slope;
        public double Intercept => _inner.Intercept;
        public int AnchorCount => _inner.AnchorCount;
        public bool IsReliable => _inner.IsReliable;
        public RtCalibrationModelType ModelType => RtCalibrationModelType.Lowess;

        public double ToMinutes(double libraryRtOrIrt) => _inner.ToMinutes(libraryRtOrIrt);
        public double GetLocalSigma(double libraryRtOrIrt) => SigmaMinutes;
        public RtCalibrationModel ToRtCalibrationModel() => _inner.ToRtCalibrationModel();
    }
}