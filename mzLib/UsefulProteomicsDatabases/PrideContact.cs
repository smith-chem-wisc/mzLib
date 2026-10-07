namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// A person associated with a PRIDE Archive project — a submitter or a lab principal
    /// investigator — as returned by the PRIDE Archive REST API (v3 <c>projects/{accession}</c>).
    /// A plain data object populated by JSON deserialization.
    /// </summary>
    /// <remarks>
    /// PRIDE returns an empty string, not a missing field, for contact details the submitter chose
    /// not to publish; <see cref="Orcid"/> is empty for most historical submissions. Every member
    /// therefore defaults to the empty string rather than null.
    /// </remarks>
    public class PrideContact
    {
        /// <summary>The person's honorific (e.g. "Professor", "Dr"), if given.</summary>
        public string Title { get; set; } = string.Empty;

        /// <summary>The person's given name(s).</summary>
        public string FirstName { get; set; } = string.Empty;

        /// <summary>The person's family name.</summary>
        public string LastName { get; set; } = string.Empty;

        /// <summary>The full display name as PRIDE composes it (first name + last name).</summary>
        public string Name { get; set; } = string.Empty;

        /// <summary>The person's institution.</summary>
        public string Affiliation { get; set; } = string.Empty;

        /// <summary>The person's contact e-mail, if published; otherwise empty.</summary>
        public string Email { get; set; } = string.Empty;

        /// <summary>The person's country.</summary>
        public string Country { get; set; } = string.Empty;

        /// <summary>The person's ORCID iD, if published; otherwise empty.</summary>
        public string Orcid { get; set; } = string.Empty;

        /// <summary>PRIDE's internal identifier for this person.</summary>
        /// <remarks>
        /// PRIDE emits this value twice, as both "id" and "identifier"; they have been observed to
        /// agree. Only one is modeled here, bound to "id".
        /// </remarks>
        public string Id { get; set; } = string.Empty;
    }
}
