namespace PredictionClients.Koina.Util
{
    public enum IncompatibleParameterHandlingMode
    {
        ThrowException, // Throw exception if incompatible parameter is present.
        ReturnNull // Return null/false if parameter is not allowed. Simply skip prediction for that peptide.

        // Note: We could add more modes in the future, like UseDefault. 
    }
}