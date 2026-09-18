namespace MassSpectrometry.Deconvolution.Consensus
{
    /// <summary>
    /// One envelope inside a <see cref="CorrectedTrace"/>, carrying both its
    /// original (algorithm-reported) mass and the post-correction mass.
    /// <see cref="WasCorrected"/> is true iff <see cref="TraceCorrector"/>
    /// flagged this envelope as an off-by-one outlier and snapped its
    /// mass to the trace consensus.
    /// </summary>
    public sealed class CorrectedEnvelope
    {
        public int ScanIndex;
        public int ScanNumber;
        public double RT;
        public double OriginalMass;
        public double CorrectedMass;
        public int Charge;
        public double Intensity;
        public bool WasCorrected;

        /// <summary>
        /// Quality score of the envelope this observation came from, carried through from
        /// <see cref="MassTrace"/>. <see cref="double.NaN"/> when the trace was built without a
        /// scorer; treat NaN as missing, not as a low score.
        /// </summary>
        public double Score = double.NaN;
    }
}
