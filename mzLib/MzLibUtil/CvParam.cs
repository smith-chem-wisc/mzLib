namespace MzLibUtil
{
    /// <summary>
    /// A controlled-vocabulary parameter (cvParam), the PSI standard unit shared across proteomics
    /// data formats (mzML, mzIdentML, mzTab) and web APIs (e.g. PRIDE Archive). A cvParam is a
    /// reference to a term in a controlled vocabulary: the vocabulary (<see cref="CvLabel"/>), the
    /// term's stable identifier (<see cref="Accession"/>), its human-readable <see cref="Name"/>, an
    /// optional <see cref="Value"/>, and — when the value is numeric — an optional unit term
    /// (<see cref="UnitCvLabel"/>/<see cref="UnitAccession"/>/<see cref="UnitName"/>).
    ///
    /// Consumers should match on <see cref="Accession"/> (the stable machine identifier), NOT on
    /// <see cref="Name"/> (a display label that can change). For example, a PRIDE FTP file location is
    /// the cvParam whose <see cref="Accession"/> is "PRIDE:0000469".
    ///
    /// This type is intentionally serialization-agnostic: it carries no JSON/XML attributes and takes
    /// no serializer dependency. JSON deserializers (Newtonsoft, System.Text.Json) map wire fields
    /// such as "cvLabel"/"accession"/"name"/"value" onto these properties by case-insensitive name
    /// matching, so no annotations are required.
    /// </summary>
    /// <remarks>
    /// <see cref="CvLabel"/> corresponds to the "cvRef" attribute in the mzML/mzIdentML XML schemas and
    /// to "cvLabel" in PRIDE's JSON — both name the controlled vocabulary the term comes from.
    /// </remarks>
    public record CvParam
    {
        /// <summary>
        /// Parameterless constructor. Used by object initializers and by deserializers, which then set
        /// the individual properties.
        /// </summary>
        public CvParam()
        {
        }

        /// <summary>
        /// Convenience constructor. Unit fields are optional and default to empty (no unit).
        /// </summary>
        /// <param name="cvLabel">The controlled vocabulary the term is from (e.g. "PRIDE", "MS").</param>
        /// <param name="accession">The stable term accession (e.g. "PRIDE:0000469").</param>
        /// <param name="name">The term's human-readable name (e.g. "FTP Protocol").</param>
        /// <param name="value">The parameter value, if any.</param>
        /// <param name="unitCvLabel">The controlled vocabulary of the unit term, if any.</param>
        /// <param name="unitAccession">The unit term accession, if any.</param>
        /// <param name="unitName">The unit term name, if any.</param>
        public CvParam(string cvLabel, string accession, string name, string value,
            string unitCvLabel = "", string unitAccession = "", string unitName = "")
        {
            CvLabel = cvLabel;
            Accession = accession;
            Name = name;
            Value = value;
            UnitCvLabel = unitCvLabel;
            UnitAccession = unitAccession;
            UnitName = unitName;
        }

        /// <summary>
        /// The controlled vocabulary the term is from (the "cvRef" / "cvLabel"), e.g. "PRIDE" or "MS".
        /// </summary>
        public string CvLabel { get; init; } = string.Empty;

        /// <summary>
        /// The stable, machine-readable term identifier, e.g. "PRIDE:0000469". Match on this, not Name.
        /// </summary>
        public string Accession { get; init; } = string.Empty;

        /// <summary>
        /// The human-readable term name, e.g. "FTP Protocol".
        /// </summary>
        public string Name { get; init; } = string.Empty;

        /// <summary>
        /// The parameter value as a string. Numeric values are interpreted via the unit term, if present.
        /// </summary>
        public string Value { get; init; } = string.Empty;

        /// <summary>
        /// The controlled vocabulary of the unit term, if the value is numeric; otherwise empty.
        /// </summary>
        public string UnitCvLabel { get; init; } = string.Empty;

        /// <summary>
        /// The accession of the unit term, if the value is numeric; otherwise empty.
        /// </summary>
        public string UnitAccession { get; init; } = string.Empty;

        /// <summary>
        /// The name of the unit term, if the value is numeric; otherwise empty.
        /// </summary>
        public string UnitName { get; init; } = string.Empty;
    }
}
