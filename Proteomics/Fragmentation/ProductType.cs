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

        a,
        aStar,
        aDegree,
        b,
        bStar,
        bDegree,
        //BnoB1ions,
        c,
        x,
        y,
        yStar,
        yDegree,
        zPlusOne,//This is zDot plus H
        zDot,
        M, //this is the molecular ion // [M]
        D, //this is a diagnostic ion // Modification loss mass
        Ycore, //Glyco core Y ions // [pep] + Neutral core Glycan mass (such as: [pep] + [N]) //Which already consider the loss of H2O and H-transfer
        Y //Glyco Y ions // [pep] + other Glycan mass 
    }

}