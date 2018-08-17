using System.Collections.Generic;

namespace Proteomics.Fragmentation
{
    public enum ProductType
    {
        //Ion Type      Neutral Mr
        //a             [N]+[M]-CHO
        //a*	        a-NH3
        //a°	        a-H2O
        //b             [N]+[M]-H
        //b*	        b-NH3
        //b°	        b-H2O
        //c             [N]+[M]+NH2
        //d             a – partial side chain
        //v             y – complete side chain
        //w             z – partial side chain
        //x             [C]+[M]+CO-H
        //y             [C]+[M]+H
        //y*	        y-NH3
        //y°	        y-H2O
        //z             [C]+[M]-NH2

        A,
        Astar,
        Adot,
        B,
        Bstar,
        Bdot,
        //BnoB1ions,
        C,
        X,
        Y,
        Ystar,
        Ydot,
        Zdot,// A Zdot ion is also known as z+1. It is not a z-ion in the Biemann nomenclature. It differs from a y-ion by N-1 H-1;
        M, //this is the molecular ion
        D //this is a diagnostic ion
    }

}