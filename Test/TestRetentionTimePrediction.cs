using NUnit.Framework;
using Proteomics;
using Proteomics.RetentionTimePrediction;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;

namespace Test
{
    /*
     * Original author: Brendan MacLean <brendanx .at. u.washington.edu>,
     *                  MacCoss Lab, Department of Genome Sciences, UW
     *
     * Copyright 2009 University of Washington - Seattle, WA
     * 
     * Licensed under the Apache License, Version 2.0 (the "License");
     * you may not use this file except in compliance with the License.
     * You may obtain a copy of the License at
     *
     *     http://www.apache.org/licenses/LICENSE-2.0
     *
     * Unless required by applicable law or agreed to in writing, software
     * distributed under the License is distributed on an "AS IS" BASIS,
     * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
     * See the License for the specific language governing permissions and
     * limitations under the License.
     */

    [TestFixture]
    public sealed class TestRetentionTimePrediction
    {
        private readonly object[,] _peptides300A;
        private readonly object[,] _peptides100A;
        private readonly object[,] _peptidesCZE;

        public TestRetentionTimePrediction()
        {
            // These peptides were taken from the supporting information for the
            // Krokhin, Anal. Chem., 2006 paper

            // Sequence-Specific Retention Calculator. Algorithm for Peptide
            //    Retention Prediction in Ion-Pair RP-HPLC: Application to
            //    300- and 100-Å Pore Size C18 Sorbents

            // http://tinyurl.com/Krokhin2006-AC

            // Supporting information:
            // http://tinyurl.com/Krokhin2006

            // The version currently implemented is SSRCalc v3.0, which
            // does not match v3.3 reported in the supporting information.
            // Reported values are noted in comments, where they differ.

            //original peptide values were not reproduced; values found at implementation were used below
            _peptides300A = new object[,]
            {
                {"LVEYR", 10.69},
                {"EVQPGR", 3.92},
                {"NQYLR", 10.39},
                {"HREER", 1.95},
                {"YQYQR", 5.68},
                {"NNEVQR", 3.77},
                {"NPFLFK", 27.33},
                {"EDPEER", 2.79},
                {"YETIEK", 8.39},
                {"NEQQER", 0.99},
                {"SSESQER", 1.34},
                {"EGEEEER", 2.06},
                {"EQIEELK", 14.34},
                {"NSYNLER", 11.59},
                {"QEDEEEK", 0.85},
                {"RNPFLFK", 28.86},
                {"REDPEER", 3.49},
                {"FEAFDLAK", 29.13},
                {"GNLELLGLK", 32.08},
                {"QEGEKEEK", 0.88},
                {"LFEITPEK", 24.2},
                {"VLLEEQEK", 17.1},
                {"EQIEELKK", 13.61},
                {"EKEDEEEK", 1.2},
                {"SHKPEYSNK", 6.08},
                {"LFEITPEKK", 22.79},
                {"EEDEEEGQR", 1.89},
                {"AIVVLLVNEGK", 32.71},
                {"QEDEEEKQK", 0.66},
                {"NILEASYNTR", 20.09},
                {"AILTVLSPNDR", 29.18},
                {"QQGEETDAIVK", 12.18},
                {"VLLEEQEKDR", 17.24},
                {"HGEWRPSYEK", 16.5},
                {"LVDLVIPVNGPGK", 31.14},
                {"RQQGEETDAIVK", 13.14},
                {"QSHFANAEPEQK", 11.27},
                {"SDLFENLQNYR", 30.44},
                {"SLPSEFEPINLR", 33.12},
                {"RHGEWRPSYEK", 16.4},
                {"ELTFPGSVQEINR", 28.46},
                {"KSLPSEFEPINLR", 32.53},
                {"RSDLFENLQNYR", 29.38},
                {"EEDEEQVDEEWR", 20.02},
                {"WEREEDEEQVDEEWR", 27.02},
                {"NFLSGSDDNVISQIENPVK", 34.63},
                {"LPAGTTSYLVNQDDEEDLR", 31.49},
                {"HGEWRPSYEKQEDEEEK", 17.96},
                {"HGEWRPSYEKEEDEEEGQR", 19.54},
                {"AKPHTIFLPQHIDADLILVVLSGK", 51.49},
                {"LPAGTTSYLVNQDDEEDLRLVDLVIPVNGPGK", 48.93},
                {"LSPGDVVIIPAGHPVAITASSNLNLLGFGINAENNER", 48.29},
                {"FDQR", 4.38},
                {"LLEYK", 14.65},
                {"ILENQK", 7.41},
                {"QVQNYK", 4.12},
                {"NSFNLER", 17.38},
                {"DSFNLER", 17.4},
                {"DDNEELR", 7.78},
                {"GQIEELSK", 14.38},
                {"VLLEEHEK", 16.5},
                {"FFEITPEK", 26.34},
                {"GDFELVGQR", 22.76},
                {"NENQQEQR", 0.39},
                {"GPIYSNEFGK", 21.85},
                {"AIVIVTVNEGK", 25.07},
                {"SDPQNPFIFK", 27.71},
                {"IFENLQNYR", 24.28},
                {"AILTVLKPDDR", 28.26},
                {"LPAGTIAYLVNR", 29.86},
                {"QQSQEENVIVK", 14.4},
                {"SVSSESEPFNLR", 23.84},
                {"SRGPIYSNEFGK", 21.2},
                {"EGSLLLPHYNSR", 26.13},
                {"QSHFADAQPQQR", 11.06},
                {"ELAFPGSAQEVDR", 24.71},
                {"RQQSQEENVIVK", 15.42},
                {"KSVSSESEPFNLR", 23.77},
                {"FQTLFENENGHIR", 28.5},
                {"VLLEEHEKETQHR", 16.28},
                {"NILEASFNTDYEEIEK", 35.62},
                {"KEDDEEEEQGEEEINK", 11.09},
                {"NPQLQDLDIFVNSVEIK", 42.27},
                {"ASSNLDLLGFGINAENNQR", 37.00},
                {"AILTVLKPDDRNSFNLER", 37.94},
                {"NFLAGDEDNVISQVQRPVK", 33.85},
                {"SKPHTIFLPQHTDADYILVVLSGK", 45.74},
                {"FFEITPEKNPQLQDLDIFVNSVEIK", 51.59},
                {"QVQLYR", 12.93},
                {"NPIYSNK", 9.96},
                {"DDNEDLR", 7.55},
                {"EQIEELSK", 14.5},
                {"SRNPIYSNK", 10.29},
                {"AIVIVTVTEGK", 26.18},
                {"SDQENPFIFK", 26.95},
                {"LPAGTIAYLANR", 27.05},
                {"SVSSESGPFNLR", 22.76},
                {"QEINEENVIVK", 21.36},
                {"EGSLLLPNYNSR", 26.4},
                {"QSYFANAQPLQR", 23.73},
                {"ELAFPGSSHEVDR", 22.94},
                {"RQEINEENVIVK", 22.8},
                {"FQTLYENENGHIR", 24.55},
                {"VLLEQQEQEPQHR", 19.09},
                {"NILEAAFNTNYEEIEK", 37.13},
                {"NQQLQDLDIFVNSVDIK", 41.34},
                {"LPAGTIAYLANRDDNEDLR", 33.2},
                {"NFLAGEEDNVISQVERPVK", 34.14},
                {"SKPHTLFLPQYTDADFILVVLSGK", 52.8},
                {"VLDLAIPVNKPGQLQSFLLSGTQNQPSLLSGFSK", 51.34},
                {"LSPGDVFVIPAGHPVAINASSDLNLIGFGINAENNER", 48.61},
                {"SFLPSK", 17.38},
                {"EGLTFR", 17.83},
                {"TILFLK", 30.69},
                {"NLFEGGIK", 24.01},
                {"DKPWWPK", 24.74},
                {"DENFGHLK", 15.61},
                {"FTPPHVIR", 23.05},
                {"DSSSPYGLR", 14.92},
                {"SSDFLAYGIK", 28.65},
                {"NNDPSLYHR", 14.24},
                {"QLSVVHPINK", 21.28},
                {"ENPHWTSDSK", 10.92},
                {"NDSELQHWWK", 27.18},
                {"SYLPSETPSPLVK", 28.38},
                {"EIFRTDGEQVLK", 26.5},
                {"SNLDPAEYGDHTSK", 14.78},
                {"SLTLEDVPNHGTIR", 26.63},
                {"LPLDVISTLSPLPVVK", 44.43},
                {"DPNSEKPATETYVPR", 16.41},
                {"VGPVQLPYTLLHPSSK", 33.89},
                {"FQTLIDLSVIEILSR", 56.36},
                {"YWVFTDQALPNDLIK", 40.64},
                {"KDPNSEKPATETYVPR", 15.78},
                {"LFILDYHDTFIPFLR", 53.07},
                {"VILPADEGVESTIWLLAK", 44.06},
                {"SLSDR", 4.42},
                {"ATLQR", 5.84},
                {"YRDR", 2.75},
                {"HIVDR", 8.12},
                {"FLVPAR", 20.89},
                {"SNNPFK", 9.3},
                {"FSYVAFK", 25.59},
                {"LDALEPDNR", 18.08},
                {"LSAEHGSLHK", 10.95},
                {"GEEEEEDKK", 1.31},
                {"GGLSIISPPEK", 24.34},
                {"QEEDEDEEK", 1.39},
                {"TVTSLDLPVLR", 31.92},
                {"ALTVPQNYAVAAK", 22.3},
                {"QEEEEDEDEER", 4.3},
                {"QEEDEDEEKQPR", 3.67},
                {"EQPQQNECQLER", 10.01},
                {"QEQENEGNNIFSGFK", 24.49},
                {"IESEGGLIETWNPNNK", 30.54},
                {"QEEEEDEDEERQPR", 5.81},
                {"LNIGPSSSPDIYNPEAGR", 26.82},
                {"LAGTSSVINNLPLDVVAATFNLQR", 44.9},
                {"FYLAGNHEQEFLQYQHQQGGK", 32.37},
                {"RFYLAGNHEQEFLQYQHQQGGK", 32.44},
                {"IEKEDVR", 7.69},
                {"VDEVFER", 18.12},
                {"GIIGLVAEDR", 28.64},
                {"QYDEEDKR", 3.82},
                {"EVAFDIAAEK", 27.09},
                {"SLWPFGGPFK", 35.79},
                {"FNLEEGDIMR", 28.00},
                {"GELETVLDEQK", 23.2},
                {"KSLWPFGGPFK", 35.46},
                {"KPESVLNTFSSK", 23.26},
                {"KSSISYHNINAK", 15.73},
                {"FGSLFEVGPSQEK", 29.86},
                {"NIENYGLAVLEIK", 35.3},
                {"EEFFFPYDNEER", 32.62},
                {"SPFNIFSNNPAFSNK", 32.81},
                {"KEEFFFPYDNEER", 32.72},
                {"EVAFDIAAEKVDEVFER", 44.39},
                {"ANAFLSPHHYDSEAILFNIK", 42.2},
                {"LYIAAFHMPPSSGSAPVNLEPFFESAGR", 44.37},
                {"EHEEEEEQEQEEDENPYVFEDNDFETK", 29.16},
                {"HKEHEEEEEQEQEEDENPYVFEDNDFETK", 26.5},
                {"QHEPR", 2.44},
                {"SPQDER", 1.8},
                {"RQQQQR", 1.77},
                {"IVNSEGNK", 5.04},
                {"HSQVAQIK", 10.92},
                {"LRSPQDER", 6.02},
                {"GDLYNSGAGR", 12.19},
                {"LSAEYVLLYR", 32.5},
                {"AAVSHVNQVFR", 23.14},
                {"ATPGEVLANAFGLR", 33.49},
                {"ISTVNSLTLPILR", 37.05},
                {"KEEEEEEQEQR", 4.03},
                {"HSEKEEEDEDEPR", 5.94},
                {"KEDEDEDEEEEEER", 6.39},
                {"GVLGLAVPGCPETYEEPR", 33.41},
                {"VFYLGGNPEIEFPETQQK", 37.06},
                {"VESEAGLTETWNPNHPELK", 31.39},
                {"VEDGLHIISPELQEEEEQSHSQR", 28.77},
                {"TIDPNGLHLPSYSPSPQLIFIIQGK", 45.07},
                {"GGQQQEEESEEQNEGNSVLSGFNVEFLAHSLNTK", 37.57},
                {"RGGQQQEEESEEQNEGNSVLSGFNVEFLAHSLNTK", 36.99},
                {"ALEAFK", 16.38},
                {"TFLWGR", 26.93},
                {"NEPWWPK", 25.98},
                {"LLYPHYR", 22.29},
                {"SDYVYLPR", 25.01},
                {"EEELNNLR", 15.37},
                {"GSAEFEELVK", 26.15},
                {"SSDFLTYGLK", 29.89},
                {"ELVEVGHGDKK", 14.09},
                {"DNPNWTSDKR", 11.67},
                {"HASDELYLGER", 21.11},
                {"LPTNILSQISPLPVLK", 43.3},
                {"NWVFTEQALPADLIK", 40.97},
                {"FQTLIDLSVIEILSR", 56.36},
                {"EHLEPNLEGLTVEEAIQNK", 36.57},
                {"ATFLEGIISSLPTLGAGQSAFK", 52.05},
                {"IFFANQTYLPSETPAPLVHYR", 43.17},
                {"IYDYDVYNDLGNPDSGENHARPVLGGSETYPYPR", 36.67},
                {"SQIVR", 8.97},
                {"VEGGLR", 8.67},
                {"SEFDR", 7.5},
                {"HSYPVGR", 10.87},
                {"EQSHSHSHR", -0.82},
                {"TANSLTLPVLR", 29.66},
                {"AAVSHVQQVLR", 23.22},
                {"ENIADAAGADLYNPR", 27.31},
                {"EEEEEEEEDEEKQR", 5.84},
                {"IRENIADAAGADLYNPR", 28.95},
                {"VESEAGLTETWNPNNPELK", 31.91},
                {"VFYLGGNPETEFPETQEEQQGR", 32.3},
                {"TIDPNGLHLPSFSPSPQLIFIIQGK", 48.01},
                {"GQLVVVPQNFVVAEQAGEEEGLEYVVFK", 48.85},
                {"KGQLVVVPQNFVVAEQAGEEEGLEYVVFK", 47.37},
                {"LLENQK", 8.32},
                {"QIEELSK", 12.03},
                {"NQVQSYK", 6.05},
                {"FFEITPK", 25.11},
                {"NENQQGLR", 6.3},
                {"KQIEELSK", 13.2},
                {"ILLEEHEK", 18.62},
                {"EEDDEEEEQR", 4.04},
                {"DLTFPGSAQEVDR", 24.13},
                {"QSYFANAQPQQR", 15.52},
                {"ILLEEHEKETHHR", 17.28},
                {"NFLAGEEDNVISQIQK", 32.48},
                {"LTPGDVFVIPAGHPVAVR", 37.28},
                {"EEDDEEEEQREEETK", 5.89},
                {"ASSNLNLLGFGINAENNQR", 35.42},
                {"NPQLQDLDIFVNYVEIK", 46.41},
                {"KNPQLQDLDIFVNYVEIK", 45.53},
                {"NENQQGLREEDDEEEEQR", 10.37},
                {"GDQYAR", 3.5},
                {"GDYYAR", 7.6},
                {"EVYLFK", 24.15},
                {"GKEVYLFK", 25.17},
                {"VLYGPTPVR", 23.15},
                {"TGYINAAFR", 23.93},
                {"TNEVYFFK", 28.18},
                {"TLDYWPSLR", 32.85},
                {"KTLDYWPSLR", 32.13},
                {"VLYGPTPVRDGFK", 27.02},
                {"YVLLDYAPGTSNDK", 31.2},
                {"SSQNNEAYLFINDK", 26.36},
                {"NTIFESGTDAAFASHK", 26.97},
            };

            //original peptide values were not reproduced; values found at implementation were used below
            _peptides100A = new object[,]
            {
                {"RQQQQR", 1.51},
                {"HSQVAQIK", 7.52},
                {"GDLYNSGAGR", 13.29},
                {"LSAEYVLLYR", 33.42},
                {"AAVSHVNQVFR", 22.00},
                {"ISTVNSLTLPILR", 36.62},
                {"KEEEEEEQEQR", 3.21},
                {"HSEKEEEDEDEPR", 3.62},
                {"ESHGQGEEEEELEK", 9.77},
                {"KEDEDEDEEEEEER", 5.68},
                {"VFYLGGNPEIEFPETQQK", 36.98},
                {"VESEAGLTETWNPNHPELK", 29.89},
                {"NGIYAPHWNINANSLLYVIR", 47.09},
                {"VEDGLHIISPELQEEEEQSHSQR", 26.16},
                {"TIDPNGLHLPSYSPSPQLIFIIQGK", 44.28},
                {"GGQQQEEESEEQNEGNSVLSGFNVEFLAHSLNTK", 36.85},
                {"SQIVR", 8.26},
                {"VEGGLR", 8.49},
                {"SEFDR", 8.04},
                {"HSYPVGR", 8.69},
                {"LSAEYVR", 14.85},
                {"QQQGDSHQK", -0.52},
                {"EQSHSHSHR", -1.04},
                {"AAVSHVQQVLR", 20.59},
                {"IVNFQGDAVFDNK", 28.05},
                {"EEEEEEEEDEEK", 6.56},
                {"ENIADAAGADLYNPR", 26.97},
                {"IRENIADAAGADLYNPR", 28.75},
                {"LNQCQLDNINALEPDHR", 28.48},
                {"VESEAGLTETWNPNNPELK", 31.4},
                {"VFYLGGNPETEFPETQEEQQGR", 32.48},
                {"TIDPNGLHLPSFSPSPQLIFIIQGK", 47.36},
                {"GQLVVVPQNFVVAEQAGEEEGLEYVVFK", 48.62},
                {"KGQLVVVPQNFVVAEQAGEEEGLEYVVFK", 46.53},
                {"ATLQR", 5.55},
                {"HIVDR", 6.17},
                {"FLVPAR", 20.14},
                {"SNNPFK", 9.39},
                {"FSYVAFK", 25.54},
                {"LDALEPDNR", 16.17},
                {"LSAEHGSLHK", 9.1},
                {"GGLSIISPPEK", 23.81},
                {"TVTSLDLPVLR", 31.9},
                {"ALTVPQNYAVAAK", 20.94},
                {"QEEEEDEDEER", 4.63},
                {"QEEDEDEEKQPR", 2.15},
                {"EQPQQNECQLER", 9.31},
                {"QEQENEGNNIFSGFK", 25.5},
                {"IESEGGLIETWNPNNK", 30.52},
                {"LNIGPSSSPDIYNPEAGR", 26.95},
                {"FYLAGNHEQEFLQYQHQQGGK", 31.12},
                {"RFYLAGNHEQEFLQYQHQQGGK", 31.15},
                {"TFLWGR", 27.87},
                {"DEAFGHLK", 17.11},
                {"NEPWWPK", 27.89},
                {"LLYPHYR", 22.39},
                {"SDYVYLPR", 24.23},
                {"EEELNNLR", 16.27},
                {"DNPNWTSDK", 12.87},
                {"ELVEVGHGDK", 12.46},
                {"KNEPWWPK", 26.5},
                {"GSAEFEELVK", 26.45},
                {"SSDFLTYGLK", 31.14},
                {"DNPNWTSDKR", 12.47},
                {"HASDELYLGER", 20.74},
                {"QDSELQAWWK", 31.32},
                {"LDSQIYGDHTSK", 12.76},
                {"LPTNILSQISPLPVLK", 42.46},
                {"NWVFTEQALPADLIK", 40.85},
                {"TWVQDYVSLYYTSDEK", 40.68},
                {"EHLEPNLEGLTVEEAIQNK", 35.5},
                {"ATFLEGIISSLPTLGAGQSAFK", 52.44},
                {"LVVEDYPYAVDGLEIWAIIK", 56.81},
                {"IFFANQTYLPSETPAPLVHYR", 42.73},
                {"LLEYK", 15.04},
                {"ILENQK", 6.86},
                {"QVQNYK", 4.56},
                {"NSFNLER", 17.52},
                {"DDNEELR", 7.78},
                {"GQIEELSK", 13.32},
                {"VLLEEHEK", 15.51},
                {"FFEITPEK", 26.44},
                {"GDFELVGQR", 21.8},
                {"NENQQEQR", 0.98},
                {"GPIYSNEFGK", 22.5},
                {"AIVIVTVNEGK", 23.63},
                {"SDPQNPFIFK", 28.93},
                {"IFENLQNYR", 25.32},
                {"AILTVLKPDDR", 26.48},
                {"LPAGTIAYLVNR", 29.64},
                {"QQSQEENVIVK", 13.52},
                {"SVSSESEPFNLR", 24.63},
                {"EGSLLLPHYNSR", 26.35},
                {"QSHFADAQPQQR", 10.22},
                {"ELAFPGSAQEVDR", 24.87},
                {"RQQSQEENVIVK", 12.27},
                {"KSVSSESEPFNLR", 23.17},
                {"VLLEEHEKETQHR", 12.63},
                {"EDDEEEEQGEEEINK", 12.29},
                {"NILEASFNTDYEEIEK", 36.78},
                {"KEDDEEEEQGEEEINK", 9.68},
                {"NPQLQDLDIFVNSVEIK", 41.82},
                {"ASSNLDLLGFGINAENNQR", 37.07},
                {"AILTVLKPDDRNSFNLER", 36.08},
                {"NFLAGDEDNVISQVQRPVK", 32.62},
                {"SKPHTIFLPQHTDADYILVVLSGK", 44.02},
                {"FFEITPEKNPQLQDLDIFVNSVEIK", 50.09},
                {"ATLTVLK", 19.31},
                {"QVQLYR", 13.55},
                {"NPIYSNK", 10.76},
                {"DDNEDLR", 7.65},
                {"EQIEELSK", 13.66},
                {"SRNPIYSNK", 9.94},
                {"AIVIVTVTEGK", 24.8},
                {"SDQENPFIFK", 28.48},
                {"LPAGTIAYLANR", 27.04},
                {"SVSSESGPFNLR", 23.5},
                {"QEINEENVIVK", 20.6},
                {"EGSLLLPNYNSR", 26.9},
                {"KSVSSESGPFNLR", 22.18},
                {"QSYFANAQPLQR", 24.51},
                {"ELAFPGSSHEVDR", 22.14},
                {"RQEINEENVIVK", 20.66},
                {"VLLEQQEQEPQHR", 17.45},
                {"NILEAAFNTNYEEIEK", 38.16},
                {"NQQLQDLDIFVNSVDIK", 41.05},
                {"NFLAGEEDNVISQVERPVK", 33.19},
                {"SKPHTLFLPQYTDADFILVVLSGK", 52.32},
                {"VLDLAIPVNKPGQLQSFLLSGTQNQPSLLSGFSK", 50.47},
                {"LSPGDVFVIPAGHPVAINASSDLNLIGFGINAENNER", 47.78},
                {"ETHHR", 1.61},
                {"LLENQK", 7.77},
                {"SEPFNLK", 19.08},
                {"QIEELSK", 12.13},
                {"FFEITPK", 25.15},
                {"NENQQGLR", 6.76},
                {"KQIEELSK", 11.93},
                {"ILLEEHEK", 17.73},
                {"EEDDEEEEQR", 4.18},
                {"DLTFPGSAQEVDR", 24.64},
                {"QSYFANAQPQQR", 16.28},
                {"LTPGDVFVIPAGHPVAVR", 35.56},
                {"EEDDEEEEQREEETK", 4.69},
                {"ASSNLNLLGFGINAENNQR", 35.43},
                {"NPQLQDLDIFVNYVEIK", 46.24},
                {"KNPQLQDLDIFVNYVEIK", 44.2},
                {"NENQQGLREEDDEEEEQR", 10.16},
                {"IILGPK", 17.14},
                {"GDQYAR", 3.71},
                {"GDYYAR", 8.66},
                {"EVYLFK", 24.63},
                {"GKEVYLFK", 24.92},
                {"VLYGPTPVR", 22.32},
                {"TGYINAAFR", 23.93},
                {"TNEVYFFK", 28.54},
                {"TLDYWPSLR", 32.73},
                {"KTLDYWPSLR", 32.47},
                {"VLYGPTPVRDGFK", 26.17},
                {"SSQNNEAYLFINDK", 27.06},
                {"NTIFESGTDAAFASHK", 28.26},
                {"IADMFPFFEGTVFENGIDAAYR", 53.44},
                {"VLILNK", 20.07},
                {"KEEHGR", 1.34},
                {"VDEVFER", 17.16},
                {"GIIGLVAEDR", 28.31},
                {"QYDEEDKR", 4.74},
                {"EVAFDIAAEK", 27.31},
                {"SLWPFGGPFK", 35.73},
                {"SSISYHNINAK", 15.96},
                {"GELETVLDEQK", 22.83},
                {"KSLWPFGGPFK", 36.4},
                {"KPESVLNTFSSK", 22.46},
                {"KSSISYHNINAK", 12.87},
                {"FGSLFEVGPSQEK", 29.91},
                {"NIENYGLAVLEIK", 35.07},
                {"GSMSTIHYNTNANK", 13.54},
                {"SPFNIFSNNPAFSNK", 33.56},
                {"KEEFFFPYDNEER", 32.62},
                {"ANAFLSPHHYDSEAILFNIK", 40.41},
                {"LYIAAFHMPPSSGSAPVNLEPFFESAGR", 44.56},
                {"EHEEEEEQEQEEDENPYVFEDNDFETK", 29.35},
                {"HKEHEEEEEQEQEEDENPYVFEDNDFETK", 25.79},
                {"TILFLK", 29.49},
                {"IFFANK", 21.36},
                {"NLFEGGIK", 23.81},
                {"DKPWWPK", 26.4},
                {"DENFGHLK", 15.92},
                {"FTPPHVIR", 21.02},
                {"SSDFLAYGIK", 28.9},
                {"NNDPSLYHR", 14.69},
                {"QLSVVHPINK", 19.79},
                {"ENPHWTSDSK", 10.31},
                {"HASDEVYLGQR", 16.48},
                {"NYMQVEFFLK", 38.88},
                {"NDSELQHWWK", 28.39},
                {"SYLPSETPSPLVK", 28.32},
                {"SNLDPAEYGDHTSK", 13.92},
                {"SLTLEDVPNHGTIR", 25.42},
                {"LPLDVISTLSPLPVVK", 43.49},
                {"DPNSEKPATETYVPR", 14.2},
                {"VGPVQLPYTLLHPSSK", 33.2},
                {"FQTLIDLSVIEILSR", 56.59},
                {"YWVFTDQALPNDLIK", 40.25},
                {"KDPNSEKPATETYVPR", 12.71},
                {"LFILDYHDTFIPFLR", 53.44},
                {"VILPADEGVESTIWLLAK", 43.72},
                {"TWVQDYVSLYYATDNDIK", 43.23},
                {"LAGTSSVINNLPLDVVAATFNLQR", 44.62},
                {"SSNNQLDQMPR", 13.76},
                {"LLIEDYPYAVDGLEIWTAIK", 60.37},
                {"ATPAEVLANAFGLR", 35.25},
                {"ATPGEVLANAFGLR", 33.21},
                {"LVEYR", 10.95},
                {"NQYLR", 10.55},
                {"HREER", 1.79},
                {"YQYQR", 6.35},
                {"NNEVQR", 3.67},
                {"NPFLFK", 28.14},
                {"YETIEK", 9.24},
                {"NEQQER", 1.45},
                {"SSESQER", 1.42},
                {"EGEEEER", 2.81},
                {"EQIEELK", 13.84},
                {"NSYNLER", 12.12},
                {"QEDEEEK", 1.89},
                {"RNPFLFK", 28.6},
                {"REDPEER", 3.91},
                {"FEAFDLAK", 29.75},
                {"GNLELLGLK", 30.94},
                {"QEGEKEEK", 0.99},
                {"LFEITPEK", 23.49},
                {"VLLEEQEK", 16.41},
                {"EQIEELKK", 11.68},
                {"SHKPEYSNK", 4.59},
                {"AIVVLLVNEGK", 31.28},
                {"NILEASYNTR", 21.03},
                {"AILTVLSPNDR", 28.31},
                {"QQGEETDAIVK", 11.83},
                {"VLLEEQEKDR", 15.36},
                {"HGEWRPSYEK", 15.26},
                {"RQQGEETDAIVK", 10.42},
                {"QSHFANAEPEQK", 10.35},
                {"SDLFENLQNYR", 32.1},
                {"SLPSEFEPINLR", 33.7},
                {"RHGEWRPSYEK", 13.76},
                {"ELTFPGSVQEINR", 28.46},
                {"KSLPSEFEPINLR", 32.00},
                {"RSDLFENLQNYR", 29.06},
                {"EEDEEQVDEEWR", 20.46},
                {"WEREEDEEQVDEEWR", 26.97},
                {"NFLSGSDDNVISQIENPVK", 34.52},
                {"LPAGTTSYLVNQDDEEDLR", 31.39},
                {"HGEWRPSYEKEEDEEEGQR", 17.62},
                {"AKPHTIFLPQHIDADLILVVLSGK", 49.06},
                {"LSPGDVVIIPAGHPVAITASSNLNLLGFGINAENNER", 47.19},
            };

            _peptidesCZE = new object[,]
            {
                {"DDDRDD", 13.69},
                {"EEEKEE", 15.42},
                {"NNNHNN", 17.15},
                {"QQQGQQ", 10.88},
                {"KKKKKK", 33.92},
                {"EDNHKRQM", 21.86},
                {"QNEHKRDE", 22.20}
            };
        }

        /// <summary>
        ///A test for ScoreSequence with 300A column
        ///</summary>
        [Test]
        public void SSRCalc3_300A_Test()
        {
            SSRCalc3 calc = new SSRCalc3("SSRCalc 3.0 (300A)", SSRCalc3.Column.A300);

            for (int i = 0; i < _peptides300A.GetLength(0); i++)
            {
                var peptide = new PeptideWithSetModifications((string)_peptides300A[i, 0], new Dictionary<string, Modification>());
                double expected = (double)_peptides300A[i, 1];
                double actual = calc.ScoreSequence(peptide);

                // Round the returned value to match the values presented
                // in the supporting information of the SSRCalc 3 publication.
                // First cast to float, since the added precision of the double
                // caused the double representation of 12.805 to round to 12.80
                // instead of 12.81.  When diffed with 12.81 the result was
                // 0.005000000000002558.
                double actualRound = Math.Round((float)actual, 2);

                // Extra conditional added to improve debugging of issues.
                if (Math.Abs(expected - actual) > 0.005)
                {
                    Assert.AreEqual(expected, actualRound, "Peptide {0}", peptide);
                }
            }
        }

        /// <summary>
        ///A test for ScoreSequence with 100A column
        ///</summary>
        // Problems with the results from the article
        [Test]
        public void SSRCalc3_100A_Test()
        {
            SSRCalc3 calc = new SSRCalc3("SSRCalc 3.0 (100A)", SSRCalc3.Column.A100);

            for (int i = 0; i < _peptides100A.GetLength(0); i++)
            {
                var peptide = new PeptideWithSetModifications((string)_peptides100A[i, 0], new Dictionary<string, Modification>());
                object obj = _peptides100A[i, 1];
                double expected = (double)_peptides100A[i, 1];
                double actual = calc.ScoreSequence(peptide);

                // Round the returned value to match the values presented
                // in the supporting information of the SSRCalc 3 publication.
                // First cast to float, since the added precision of the double
                // caused the double representation of 12.805 to round to 12.80
                // instead of 12.81.  When diffed with 12.81 the result was
                // 0.005000000000002558.
                double actualRound = Math.Round((float)actual, 2);

                // Extra conditional added to improve debugging of issues.
                if (Math.Abs(expected - actual) > 0.005)
                {
                    Assert.AreEqual(expected, actualRound, "Peptide {0}", peptide);
                }
            }
        }

        /// <summary>
        ///A test for CZE retention time prediction
        ///</summary>
        [Test]
        public void CZE_RetentionTime_Test()
        {
            CZE testCZE = new CZE(1,1);

            double expElutionTime = 1;
            double expMu = Math.Round(testCZE.ExperimentalElectrophoreticMobility(expElutionTime), 0);
            Assert.AreEqual(expMu, 16666667);
            double theoreticalElutionTime = testCZE.TheoreticalElutionTime(expMu);
            Assert.AreEqual(Math.Round(expElutionTime, 5), Math.Round(theoreticalElutionTime, 5));


            for (int i = 0; i < _peptidesCZE.GetLength(0); i++)
            {
                var peptide = new PeptideWithSetModifications((string)_peptidesCZE[i, 0], new Dictionary<string, Modification>());
                object obj = _peptidesCZE[i, 1];
                double expected = (double)_peptidesCZE[i, 1];
                double actual = CZE.PredictedElectrophoreticMobility(peptide.BaseSequence, peptide.MonoisotopicMass);

                // Round the returned value to match the values presented
                // in the supporting information of the SSRCalc 3 publication.
                // First cast to float, since the added precision of the double
                // caused the double representation of 12.805 to round to 12.80
                // instead of 12.81.  When diffed with 12.81 the result was
                // 0.005000000000002558.
                double actualRound = Math.Round((float)actual, 2);

                // Extra conditional added to improve debugging of issues.
                if (Math.Abs(expected - actual) > 0.005)
                {
                    Assert.AreEqual(expected, actualRound, "Peptide {0}", peptide);
                }
            }
        }
    }
}



