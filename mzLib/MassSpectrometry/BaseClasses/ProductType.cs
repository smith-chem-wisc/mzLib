using Chemistry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
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

        // Base ions are for Nucleic acids in which the base is cleaved as a neutral loss during fragmentation
        // schematic for RNA fragmentation modes can be found below
        // https://www.researchgate.net/figure/The-standard-nomenclature-for-oligonucleotide-fragmentation-during-collisioninduced_fig6_271536997


        a,
        aStar,
        aDegree,
        aWaterLoss,
        aBaseLoss,          
        b,
        bAmmoniaLoss,
        bWaterLoss,
        //BnoB1ions,            
        bBaseLoss,
        c,
        cWaterLoss,
        cBaseloss,
        d,
        dWaterLoss,
        dBaseLoss,
        w,
        wWaterLoss,
        wBaseLoss,
        x,
        xWaterLoss,
        xBaseLoss,
        y,
        yAmmoniaLoss,
        yWaterLoss,
        yBaseLoss,
        z,
        zPlusOne,       //This is zDot plus H
        zDot,
        zWaterLoss,
        zBaseLoss,
        M,              //this is the molecular ion // [M]
        D,              //this is a diagnostic ion // Modification loss mass
        Ycore,          //Glyco core Y ions // [pep] + Neutral core Glycan mass (such as: [pep] + [N]) //Which already consider the loss of H2O and H-transfer
        Y               //Glyco Y ions // [pep] + other Glycan mass 
    }
}
