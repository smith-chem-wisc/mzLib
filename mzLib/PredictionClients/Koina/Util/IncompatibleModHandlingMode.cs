namespace PredictionClients.Koina.Util
{
    /// <summary>
    /// Mode for handling incompatible modifications during model inference.
    /// </summary>
    public enum IncompatibleModHandlingMode
    {
        RemoveIncompatibleMods, // Strip unsupported modifications and predict with remaining mods
        UsePrimarySequence, // Ignore all modifications and use only base sequence
        ThrowException, // Throw exception if incompatible modifications present.
        ReturnNull // Return null if incompatible modifications present. Simply skip prediction for that peptide.
    }
}
