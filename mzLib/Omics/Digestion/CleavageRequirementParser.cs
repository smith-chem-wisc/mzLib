using MzLibUtil;
using Omics.Modifications;
using System;
using System.Collections.Generic;

namespace Omics.Digestion
{
    /// <summary>
    /// Reads the compact text form of a <see cref="CleavageRequirement"/> used by the protease data
    /// files -- <c>P2:O-glycan</c>, <c>P1':O-glycan</c>, <c>P4:O-glycan</c> -- so a glycoprotease can be
    /// configured from <c>proteases.tsv</c> rather than only in code.
    /// </summary>
    /// <remarks>
    /// <para><b>The syntax deliberately reads as the enzymology does.</b> A subsite is written the way
    /// every paper writes it: <c>P</c> then the subsite number, with a trailing apostrophe for the prime
    /// side. <c>P2</c> is two residues before the severed bond, <c>P1'</c> is the residue immediately
    /// after it. So StcE's rule is <c>P2:O-glycan</c> and the OgpA family's is <c>P1':O-glycan</c>,
    /// which is exactly how the family table in the literature prints them.</para>
    ///
    /// <para><b>Why the class names are coarse.</b> The value after the colon names a glycosylation
    /// CLASS, not a structure, because that is all a <see cref="Modification"/> can currently be
    /// resolved to. Real enzymes discriminate far more finely -- OpeRATOR needs core 1 and is refused by
    /// the Tn antigen, AMUC_1438 needs Tn and is refused by anything larger -- and none of that can be
    /// written here yet. A rule written at class level admits glycoforms the real enzyme would refuse,
    /// which costs search space and never an identification.</para>
    ///
    /// <para>Whitespace is ignored and the keyword is case-insensitive, because these values are typed
    /// by hand into a tab-separated file.</para>
    /// </remarks>
    public static class CleavageRequirementParser
    {
        /// <summary>
        /// Parses one requirement, or returns null when <paramref name="text"/> is empty -- which is the
        /// case for every protease that cleaves on sequence alone, i.e. nearly all of them.
        /// </summary>
        /// <exception cref="MzLibException">
        /// When the text is present but malformed. A silently ignored requirement would let a
        /// glycoprotease digest as if it had none, which is the over-digestion this feature exists to
        /// stop, so a typo has to be loud.
        /// </exception>
        /// <summary>
        /// Every condition in a "Cleavage Requirement" cell, which may name more than one, separated by
        /// semicolons. An empty cell yields an empty list, which is the common case.
        /// </summary>
        /// <remarks>
        /// IMPa is why this is a list: <c>P1':O-glycan;P1:!O-glycan</c> says it cleaves N-terminal to a
        /// glycosylated Ser/Thr but not when the residue before the bond is itself glycosylated. Those are
        /// two conditions on two residues, and both must hold. The <c>!</c> marks the forbidding one.
        /// </remarks>
        public static List<CleavageRequirement> ParseAll(string text)
        {
            var requirements = new List<CleavageRequirement>();
            if (string.IsNullOrWhiteSpace(text))
            {
                return requirements;
            }

            foreach (string clause in text.Split(';'))
            {
                if (string.IsNullOrWhiteSpace(clause))
                {
                    continue;
                }

                requirements.Add(Parse(clause));
            }

            return requirements;
        }

        public static CleavageRequirement Parse(string text)
        {
            if (string.IsNullOrWhiteSpace(text))
            {
                return null;
            }

            string trimmed = text.Replace(" ", string.Empty).Trim();

            int colon = trimmed.IndexOf(':');
            if (colon < 0)
            {
                throw new MzLibException("Unrecognized cleavage requirement '" + text
                    + "'. Expected a subsite and a modification class separated by a colon, for example "
                    + "P2:O-glycan for StcE or P1':O-glycan for OpeRATOR.");
            }

            string subsiteText = trimmed.Substring(0, colon);
            string classText = trimmed.Substring(colon + 1);

            if (subsiteText.Length < 2 || (subsiteText[0] != 'P' && subsiteText[0] != 'p'))
            {
                throw new MzLibException("Unrecognized subsite '" + subsiteText + "' in cleavage requirement '"
                    + text + "'. A subsite is written P followed by its number, with a trailing apostrophe "
                    + "for the prime side: P1, P2, P1', P2'.");
            }

            bool isPrimeSide = subsiteText.EndsWith("'", StringComparison.Ordinal);
            string numberText = subsiteText.Substring(1, subsiteText.Length - 1 - (isPrimeSide ? 1 : 0));

            if (!int.TryParse(numberText, out int subsite) || subsite < 1)
            {
                throw new MzLibException("Unrecognized subsite number in cleavage requirement '" + text
                    + "'. Subsites are numbered outward from the severed bond starting at 1.");
            }

            // A leading '!' forbids the class at this subsite instead of requiring it.
            bool isForbidden = classText.StartsWith("!", StringComparison.Ordinal);
            if (isForbidden)
            {
                classText = classText.Substring(1);
            }

            GlycosylationClass requiredClass = ParseClass(classText, text);

            if (isForbidden)
            {
                return isPrimeSide
                    ? CleavageRequirement.PrimeForbidden(subsite, requiredClass)
                    : CleavageRequirement.NonPrimeForbidden(subsite, requiredClass);
            }

            return isPrimeSide
                ? CleavageRequirement.Prime(subsite, requiredClass)
                : CleavageRequirement.NonPrime(subsite, requiredClass);
        }

        private static GlycosylationClass ParseClass(string classText, string wholeText) =>
            classText.ToLowerInvariant() switch
            {
                "o-glycan" or "oglycan" or "o-linked" => GlycosylationClass.OLinked,
                "n-glycan" or "nglycan" or "n-linked" => GlycosylationClass.NLinked,
                _ => throw new MzLibException("Unrecognized modification class '" + classText
                    + "' in cleavage requirement '" + wholeText + "'. Supported values are O-glycan and "
                    + "N-glycan. Structure-level requirements (core 1, Tn) cannot be expressed yet, "
                    + "because a Modification does not carry a monosaccharide composition."),
            };
    }
}
