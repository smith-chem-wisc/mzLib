namespace MassSpectrometry
{
    public static class DissociationTypeExtensions
    {
        /// <summary>
        /// Whether this is collision-induced dissociation of the ion-trap kind: <see cref="DissociationType.CID"/>
        /// or <see cref="DissociationType.LowCID"/>. LowCID is the same activation read out at low resolution, so
        /// code that special-cases CID nearly always has to special-case LowCID alongside it.
        ///
        /// Deliberately excludes HCD and ISCID. HCD is beam-type CID and does not share the b1 behaviour that
        /// motivates the check, and ISCID happens in the source rather than the trap.
        /// </summary>
        public static bool IsCid(this DissociationType dissociationType)
        {
            return dissociationType == DissociationType.CID || dissociationType == DissociationType.LowCID;
        }
    }
}
