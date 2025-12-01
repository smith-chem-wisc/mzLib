using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IsobaricQuant
{
    internal class TmtFileInfo
    {
        public TmtFileInfo(string fullFilePathWithExtension, string plex, int fraction, int technicalReplicate, IReadOnlyList<TmtPlexAnnotation> annotations)
        {
            FullFilePathWithExtension = fullFilePathWithExtension;
            Plex = plex ?? string.Empty;
            Fraction = fraction;
            TechnicalReplicate = technicalReplicate;
            Annotations = annotations ?? Array.Empty<TmtPlexAnnotation>();
        }

        public string FullFilePathWithExtension { get; }
        public string Plex { get; }
        public int Fraction { get; }             // 1-based
        public int TechnicalReplicate { get; }   // 1-based
        public IReadOnlyList<TmtPlexAnnotation> Annotations { get; } // All tags for this file's plex
        public override bool Equals(object obj)
        {
            if (base.Equals(obj))
            {
                return ((TmtFileInfo)obj).FullFilePathWithExtension.Equals(FullFilePathWithExtension);
            }

            return false;
        }

        public override int GetHashCode()
        {
            return FullFilePathWithExtension.GetHashCode();
        }

        public override string ToString()
        {
            return Path.GetFileName(FullFilePathWithExtension);
        }
    }
}
