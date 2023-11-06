using MathNet.Numerics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Proteomics.ProteolyticDigestion
{
    public class ChronologerTensorResidueDictionary
    {
        private static readonly Dictionary<(char,string), int> _residueWithModToTensorInt =
            new Dictionary<(char, string), int>();
        public ChronologerTensorResidueDictionary()
        {
            FillTensorIntDictionary();
        }

        private void FillTensorIntDictionary()
        {
            _residueWithModToTensorInt.Clear();
            _residueWithModToTensorInt.Add(('A',""),1);//'Alanine
            _residueWithModToTensorInt.Add(('C', ""),2);//'Cysteine
            _residueWithModToTensorInt.Add(('D', ""),3);//'Aspartate
            _residueWithModToTensorInt.Add(('E', ""),4);//'Glutamate
            _residueWithModToTensorInt.Add(('F', ""),5);//'Phenylalaline
            _residueWithModToTensorInt.Add(('G', ""),6);//'Glycine
            _residueWithModToTensorInt.Add(('H', ""),7);//'Histidine
            _residueWithModToTensorInt.Add(('I', ""),8);//'Isoleucine
            _residueWithModToTensorInt.Add(('K', ""),9);//'Lysine
            _residueWithModToTensorInt.Add(('L', ""),10);//'Leucine
            _residueWithModToTensorInt.Add(('M', ""),11);//'Methionine
            _residueWithModToTensorInt.Add(('N', ""),12);//'Asparagine
            _residueWithModToTensorInt.Add(('P', ""),13);//'Proline
            _residueWithModToTensorInt.Add(('Q', ""),14);//'Glutamine
            _residueWithModToTensorInt.Add(('R', ""),15);//'Argenine
            _residueWithModToTensorInt.Add(('S', ""),16);//'Serine
            _residueWithModToTensorInt.Add(('T', ""),17);//'Threonine
            _residueWithModToTensorInt.Add(('V', ""),18);//'Valine
            _residueWithModToTensorInt.Add(('W', ""),19);//'Tryptophane
            _residueWithModToTensorInt.Add(('Y', ""),20);//'Tyrosine
            _residueWithModToTensorInt.Add(('C', "Carbamidomethyl on C"),21);//'Carbamidomethyl
            _residueWithModToTensorInt.Add(('M', "Oxidation on M"),22);//'Oxidized
            //_residueWithModToTensorInt.Add(('C',null),23);//'S - carbamidomethylcysteine
            _residueWithModToTensorInt.Add(('E', "Glu to PyroGlu"),24);//'Pyroglutamate
            _residueWithModToTensorInt.Add(('S',"Phosphorylation on S"),25);//'Phosphoserine
            _residueWithModToTensorInt.Add(('T',"Phosphorylation on T"),26);//'Phosphothreonine
            _residueWithModToTensorInt.Add(('Y',"Phosphorylation on Y"),27);//'Phosphotyrosine
            _residueWithModToTensorInt.Add(('K',"Accetylation on K"),28);//'Acetylated
            _residueWithModToTensorInt.Add(('K',"Succinylation on K"),29);//'Succinylated
            _residueWithModToTensorInt.Add(('K',"Ubiquitination on K"),30);//'Ubiquitinated
            _residueWithModToTensorInt.Add(('K',"Methylation on K"),31);//'Monomethyl
            _residueWithModToTensorInt.Add(('K',"Dimethylation on K"),32);//'Dimethyl
            _residueWithModToTensorInt.Add(('K',"Trimethylation on K"),33);//'Trimethyl
            _residueWithModToTensorInt.Add(('R',"Methylation on R"),34);//'Monomethyl
            _residueWithModToTensorInt.Add(('R',"Dimethylation on R"),35);//'Dimethyl
        }

        public int TensorInt((char, string) residueWithMod)
        {
            if (_residueWithModToTensorInt.ContainsKey(residueWithMod))
            {
                return _residueWithModToTensorInt[residueWithMod];
            }

            return -1;
        }

    }
}
