using Chemistry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Transcriptomics
{
    [Flags]
    public enum FragmentType
    {
        None = 0,
        a = 1 << 0,
        adot = 1 << 1,
        aBase = 1 << 2,
        b = 1 << 3,
        bdot = 1 << 4,
        bBase = 1 << 5,
        c = 1 << 6,
        cdot = 1 << 7,
        cBase = 1 << 8,
        d = 1 << 9,
        ddot = 1 << 10,
        dBase = 1 << 11,
        dH2O = 1 << 12, // d-H20
        w = 1 << 13,
        wdot = 1 << 14,
        wBase = 1 << 15,
        x = 1 << 16,
        xdot = 1 << 17,
        xBase = 1 << 18,
        y = 1 << 19,
        ydot = 1 << 20,
        yBase = 1 << 21,
        z = 1 << 22,
        zdot = 1 << 23,
        zBase = 1 << 24,
        Internal = 1 << 25,
        All = (1 << 25) - 1, // Handy way of setting all below the 32th bit
    }

    public static class FragmentTypeExtensions
    {
        private static readonly Dictionary<FragmentType, ChemicalFormula> FragmentIonCaps = new Dictionary<FragmentType, ChemicalFormula>
        {
            {FragmentType.a, ChemicalFormula.ParseFormula("H")},
            {FragmentType.adot, new ChemicalFormula()},
            {FragmentType.b, ChemicalFormula.ParseFormula("OH")},
            {FragmentType.bdot, ChemicalFormula.ParseFormula("OH2")},
            {FragmentType.c, ChemicalFormula.ParseFormula("O3H2P")},
            {FragmentType.cdot, ChemicalFormula.ParseFormula("O3HP")},
            {FragmentType.d, ChemicalFormula.ParseFormula("O4H2P")},
            {FragmentType.ddot, ChemicalFormula.ParseFormula("O4H3P")},

            {FragmentType.w, ChemicalFormula.ParseFormula("H")},
            {FragmentType.wdot, ChemicalFormula.ParseFormula("H2")},
            {FragmentType.x, ChemicalFormula.ParseFormula("O-1H")},
            {FragmentType.xdot, ChemicalFormula.ParseFormula("O-1")},
            {FragmentType.y, ChemicalFormula.ParseFormula("O-3P-1")},
            {FragmentType.ydot, ChemicalFormula.ParseFormula("O-3HP-1")},
            {FragmentType.z, ChemicalFormula.ParseFormula("O-4P-1")},
            {FragmentType.zdot, ChemicalFormula.ParseFormula("O-4H-1P-1")},

            //fragment - Base chemical formula is the corresponding fragment chemical formula subtracing 1 H as H is lost when base is removed
            {FragmentType.aBase, ChemicalFormula.ParseFormula("H-2")}, // "H-1" -H 
            {FragmentType.bBase, ChemicalFormula.ParseFormula("O1H-2")}, //"OH1" -H

            {FragmentType.cBase, ChemicalFormula.ParseFormula("O3H-1P")}, //"O3P" -H
            {FragmentType.dBase, ChemicalFormula.ParseFormula("O4H-1P")}, //"O4H2P" -H
            {FragmentType.wBase, new ChemicalFormula()}, //"H"-H
            {FragmentType.xBase, ChemicalFormula.ParseFormula("O-1")}, //"O-1H" -H
            {FragmentType.yBase, ChemicalFormula.ParseFormula("O-3H-1P-1")}, //"O-3P-1" -H
            {FragmentType.zBase, ChemicalFormula.ParseFormula("O-4H-2P-1")}, //"O-4H-1P-1" -1
            //d-H2O
            {FragmentType.dH2O, ChemicalFormula.ParseFormula("O3P")},

        };

        public static IEnumerable<FragmentType> GetIndividualFragmentTypes(this FragmentType fragmentType)
        {
            if (fragmentType == FragmentType.None)
                yield break;
            foreach (FragmentType site in Enum.GetValues(typeof(FragmentType)))
            {
                if (site == FragmentType.None || site == FragmentType.All || site == FragmentType.Internal)
                {
                    continue;
                }
                if ((fragmentType & site) == site)
                {
                    yield return site;
                }
            }
        }

        public static Terminus GetTerminus(this FragmentType fragmentType)
        {
            // Super handy: http://stackoverflow.com/questions/4624248/c-logical-riddle-with-bit-operations-only-one-bit-is-set
            if (fragmentType == FragmentType.None || (fragmentType & (fragmentType - 1)) != FragmentType.None)
            {
                throw new ArgumentException("Fragment Type must be a single value to determine the terminus", "fragmentType");
            }
            var returnValue = fragmentType >= FragmentType.w ? Terminus.ThreePrime : Terminus.FivePrime;
            return returnValue;
        }

        public static ChemicalFormula GetIonCap(this FragmentType fragmentType)
        {
            if (fragmentType == FragmentType.None || (fragmentType & (fragmentType - 1)) != FragmentType.None)
            {
                throw new ArgumentException("Fragment Type must be a single value to determine the ion cap", "fragmentType");
            }

            return FragmentIonCaps[fragmentType];
        }
    }
}
