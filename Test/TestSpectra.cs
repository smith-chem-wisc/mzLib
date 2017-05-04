// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work Copyright 2016 Stefan Solntsev
//
// This file (TestSpectra.cs) is part of MassSpectrometry.Tests.
//
// MassSpectrometry.Tests is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MassSpectrometry.Tests is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with MassSpectrometry.Tests. If not, see <http://www.gnu.org/licenses/>.

using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Spectra;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public sealed class SpectrumTestFixture
    {

        #region Private Fields

        private MzmlMzSpectrum _mzSpectrumA;

        #endregion Private Fields

        #region Public Methods

        [SetUp]
        public void Setup()
        {
            double[] mz = { 328.73795, 329.23935, 447.73849, 448.23987, 482.23792, 482.57089, 482.90393, 500.95358, 501.28732, 501.62131, 611.99377, 612.32806, 612.66187, 722.85217, 723.35345 };
            double[] intensities = { 81007096.0, 28604418.0, 78353512.0, 39291696.0, 122781408.0, 94147520.0, 44238040.0, 71198680.0, 54184096.0, 21975364.0, 44514172.0, 43061628.0, 23599424.0, 56022696.0, 41019144.0 };

            _mzSpectrumA = new MzmlMzSpectrum(mz, intensities, false);
        }

        [Test]
        public void SpectrumCount()
        {
            Assert.AreEqual(15, _mzSpectrumA.Size);
        }

        [Test]
        public void SpectrumFirstMZ()
        {
            Assert.AreEqual(328.73795, _mzSpectrumA.FirstX);
        }

        [Test]
        public void SpectrumLastMZ()
        {
            Assert.AreEqual(723.35345, _mzSpectrumA.LastX);
        }

        [Test]
        public void SpectrumBasePeakIntensity()
        {
            double basePeakIntensity = _mzSpectrumA.YofPeakWithHighestY;

            Assert.AreEqual(122781408.0, basePeakIntensity);
        }

        [Test]
        public void SpectrumTIC()
        {
            double tic = _mzSpectrumA.SumOfAllY;

            Assert.AreEqual(843998894.0, tic);
        }

        [Test]
        public void SpectrumGetIntensityFirst()
        {
            Assert.AreEqual(81007096.0, _mzSpectrumA.First().Intensity);
        }

        [Test]
        public void SpectrumGetIntensityRandom()
        {
            Assert.AreEqual(44238040.0, _mzSpectrumA[6].Intensity);
        }

        [Test]
        public void MyDeconvolutionTest()
        {
            // Should find
            // [[199.09600,270.13300,341.17000,454.25400,582.31300,653.35000,752.41800,866.46100,953.49300,1052.56100,1180.62000];[174.11200,287.19600,483.31700,570.34900,698.40800,755.42900,826.46600,957.50700,1028.54400,1159.58400,1216.60600,1287.64300,1402.67000];]

            double[] mz1 = { 101.823464108355 , 101.868928075899 , 101.931436262858 , 101.968439156395 , 102.004694362599 , 102.012018646681 , 102.016726025096 , 102.024462300157 , 102.052080954714 , 102.058108230157 , 102.075335556674 , 102.083163385286 , 102.089350879442 , 102.125957040925 , 102.138217587299 , 102.147502643265 , 102.15581112802 , 102.191478865605 , 102.200924140285 , 102.209804834734 , 102.219662100394 , 102.22741363438 , 102.237057275087 , 102.244221340455 , 102.252819744788 , 102.327306188005 , 102.436811864487 , 124.562589057293 , 159.336451176171 , 159.350947155082 , 159.406108169572 , 159.423167647912 , 159.544399808387 , 159.687452231856 , 159.702848487352 , 159.736799595855 , 159.756834564604 , 159.771971418372 , 159.90011587195 , 175.118894809794 , 183.076531098408 , 186.006412579218 , 186.316367126196 , 200.102928088662 , 223.529437098656 , 226.015970508654 , 256.166767403623 , 271.139618244504 , 283.136093659581 , 288.202911315819 , 289.212350251513 , 301.152459197012 , 307.651449051853 , 307.872123627661 , 313.190255781449 , 325.188531749696 , 336.649510444858 , 342.177361266182 , 344.206676242386 , 354.17578982368 , 366.540646210211 , 372.190415843642 , 373.433377370138 , 382.167250441073 , 384.220796590447 , 406.980337981262 , 427.266316048474 , 443.227670579715 , 453.212836550288 , 455.260096022492 , 456.265048834186 , 457.271497020547 , 464.915577655843 , 467.289317090984 , 483.291016220374 , 484.324564257837 , 485.326923052252 , 498.894701545631 , 501.171821466596 , 509.31575397261 , 512.286575149651 , 515.235331920897 , 526.298907168932 , 538.315554883084 , 542.295082567411 , 554.329186491957 , 566.295968038653 , 567.235734722021 , 571.355827625014 , 572.358705222885 , 573.358836214225 , 583.319557386676 , 584.322007734643 , 586.276370870412 , 589.188933471154 , 602.90719548456 , 603.249300586871 , 606.232420457599 , 609.328578445987 , 616.277798146893 , 625.367967120547 , 637.123443071489 , 637.68277423252 , 638.242532643455 , 638.804610411015 , 639.374256605461 , 648.243354271252 , 654.355713480137 , 655.36152079164 , 665.943280218509 , 678.04159370017 , 678.667270667839 , 679.869063614226 , 682.387274552876 , 682.818919028084 , 684.384912000427 , 687.304310599646 , 696.678234545761 , 698.375942560171 , 699.415685721253 , 700.417708819314 , 718.396445847276 , 720.377298477141 , 728.351429663863 , 730.76252294782 , 730.845592536445 , 731.421830586564 , 731.564654126155 , 732.216393338015 , 732.339441310585 , 736.401794440735 , 737.417123321544 , 738.513141398653 , 739.46450486515 , 753.424163074744 , 754.42929899354 , 755.424486093126 , 756.4366411175 , 757.440861500786 , 758.445142919772 , 768.398234629638 , 776.909541018044 , 777.408019585497 , 777.912052401711 , 786.418109719859 , 800.409750636609 , 803.925529067171 , 804.406734531331 , 804.924805558702 , 805.434087445175 , 812.427985278979 , 812.928294917451 , 813.427444877611 , 820.435441958272 , 822.407261304436 , 826.457651437253 , 827.473224460865 , 828.478726593865 , 829.474463014757 , 833.415660279057 , 840.475659776697 , 841.443014596772 , 849.517549511113 , 852.427975862332 , 860.45014214549 , 861.458573992123 , 867.46888150947 , 868.473590178362 , 868.965049640236 , 869.468105885239 , 869.968354488012 , 871.441756302424 , 895.968952620628 , 896.471520580026 , 896.988187786281 , 897.464449358686 , 904.486667757606 , 904.988259145793 , 905.485944249137 , 905.980760674548 , 931.999303303165 , 940.005234089941 , 940.509755191761 , 941.009393437526 , 941.499449078285 , 954.503593393736 , 955.499085671825 , 958.515225856621 , 959.514868562356 , 960.523239373288 , 966.005710151156 , 967.50254967463 , 970.497815649793 , 974.520495574501 , 975.014823714307 , 975.528500171229 , 976.907968042295 , 983.504889678877 , 984.522537916312 , 993.021820271685 , 1010.03332455097 , 1010.53155897562 , 1011.03846046977 , 1011.5257084683 , 1019.52858747004 , 1020.55935838312 , 1029.5526636998 , 1030.55523611916 , 1031.55744232433 , 1035.56529057377 , 1036.55101713975 , 1037.54370177561 , 1038.55616197848 , 1046.05964687726 , 1046.54121855563 , 1053.56593941828 , 1054.56875598045 , 1055.05643122888 , 1055.55851090267 , 1056.04704065091 , 1082.55350267076 , 1090.09238827594 , 1090.58018559577 , 1091.07750448491 , 1122.49965976612 , 1127.58137013332 , 1143.01058490819 , 1160.58911084679 , 1161.6047449061 , 1181.62506508611 , 1217.6067092795 , 1218.61990191079 , 1219.61246447524 , 1227.14243886812 , 1227.58433734104 , 1227.79747400781 , 1228.81615985216 , 1288.6466495864 , 1289.63994457927 , 1403.68612787081 , 1404.68125393469 , 1405.66527150105 , 1502.72741657861 , 1503.75904199151 , 1504.66200414736 , 1546.18874174752 , 1573.72658503746 , 1574.79324494254 , 1575.71793580784 , 1639.31322963142 , 1639.90930428426 , 1640.59766491653 , 1641.20997506575 , 1672.96965777172 , 1673.34014447484 , 1674.34503625084 , 1674.89472377116 , 1675.32893174579 , 1675.90193490378 , 1676.47152006252 , 1677.69394307574 , 1678.14145683312 , 1693.27147460325 , 1749.83020663825 , 1833.15651142061 , 1995.55920514699};
            double[] intensities1 = { 1149.80334472656, 844.826538085938, 992.102111816406, 1849.63000488281, 1610.8984375, 1138.74731445313, 1189.287109375, 904.736877441406, 1237.765625, 927.797241210938, 712.818115234375, 1319.0712890625, 1084.42102050781, 946.796752929688, 830.744934082031, 3594.71557617188, 1451.08483886719, 939.439025878906, 1159.98767089844, 1601.95935058594, 4802.05419921875, 1898.15014648438, 991.622131347656, 912.080322265625, 2710.07055664063, 1333.42980957031, 762.731811523438, 939.656433105469, 1119.10266113281, 1128.86218261719, 949.288269042969, 872.761779785156, 1300.87377929688, 2004.78112792969, 1414.6455078125, 1247.21362304688, 1846.08532714844, 1378.49584960938, 1176.70251464844, 11786.662109375, 976.286315917969, 888.003784179688, 1189.81677246094, 2102.10620117188, 1022.23864746094, 1031.68603515625, 1729.54870605469, 3692.4326171875, 1757.72741699219, 9542.0869140625, 1221.73913574219, 2613.1708984375, 1183.6826171875, 956.828430175781, 1152.12609863281, 1034.67919921875, 817.410522460938, 5239.68310546875, 957.250427246094, 1430.91040039063, 1350.10375976563, 3476.11108398438, 1046.40295410156, 823.159851074219, 2349.78076171875, 881.396484375, 1959.19921875, 2005.71484375, 969.892150878906, 8186.31640625, 1574.64733886719, 894.248107910156, 1162.17724609375, 1051.11645507813, 1074.27600097656, 29082.455078125, 7954.880859375, 1179.47277832031, 945.611145019531, 1026.74230957031, 1094.75, 1201.76843261719, 1860.970703125, 966.66259765625, 1404.51293945313, 2586.81079101563, 1752.26000976563, 880.1650390625, 30420.283203125, 9795.7490234375, 1974.94152832031, 8226.484375, 2786.30078125, 1706.71630859375, 854.429870605469, 1378.87561035156, 893.962951660156, 903.980712890625, 1596.94262695313, 1815.15869140625, 1575.81811523438, 1293.00756835938, 2148.68627929688, 1125.01745605469, 1977.63122558594, 1934.373046875, 1266.3486328125, 9765.9267578125, 3735.44848632813, 1051.90405273438, 3595.17700195313, 1542.74584960938, 1070.79907226563, 1346.79711914063, 1271.6513671875, 1141.42578125, 1315.08752441406, 1021.58056640625, 1249.45422363281, 6347.82470703125, 4705.8837890625, 2383.91577148438, 1225.45373535156, 973.904724121094, 1413.02258300781, 957.344482421875, 1759.99206542969, 2073.75439453125, 1529.18176269531, 1131.97888183594, 2510.73095703125, 1251.30712890625, 1242.21447753906, 1298.86633300781, 16051.482421875, 7099.68310546875, 2058.01196289063, 35310.234375, 14446.1083984375, 3163.90747070313, 1842.02014160156, 8903.7685546875, 9176.208984375, 3598.08740234375, 985.557067871094, 1244.15783691406, 2772.29223632813, 2116.92846679688, 1444.28527832031, 1216.28857421875, 8107.38623046875, 8803.044921875, 6288.58251953125, 1072.87951660156, 937.758666992188, 1384.03259277344, 21576.74609375, 10665.248046875, 1948.43542480469, 912.775939941406, 1240.486328125, 1471.90307617188, 859.984191894531, 1096.92553710938, 1910.80126953125, 1254.17126464844, 5433.5048828125, 2188.17700195313, 5730.34033203125, 8259.228515625, 3594.54150390625, 1121.21496582031, 1749.21044921875, 2697.98779296875, 1472.60083007813, 1107.88024902344, 6026.9619140625, 8961.7578125, 6048.22265625, 1357.50830078125, 2034.0712890625, 5366.1923828125, 7442.048828125, 4511.67041015625, 1684.14758300781, 4230.89111328125, 2854.49243164063, 18791.626953125, 10913.8466796875, 5037.36865234375, 1434.15075683594, 1796.27734375, 2295.21997070313, 1969.51989746094, 1267.80603027344, 1682.12145996094, 987.114807128906, 1914.62182617188, 1855.150390625, 1037.904296875, 3447.943359375, 5670.26708984375, 3029.24560546875, 1506.97680664063, 1444.01892089844, 1138.56628417969, 13085.6708984375, 10849.6162109375, 3364.58203125, 3635.86181640625, 4282.6416015625, 2295.60083007813, 1282.28649902344, 1295.41076660156, 1705.89807128906, 5115.326171875, 6346.48974609375, 4557.25439453125, 7274.11083984375, 1256.12829589844, 1482.93762207031, 1693.02575683594, 3228.69262695313, 1831.68579101563, 2542.2744140625, 1219.09619140625, 1263.09606933594, 2943.140625, 2654.11376953125, 3476.71948242188, 7522.21923828125, 6839.31982421875, 2918.58276367188, 1770.67346191406, 2147.80883789063, 3180.57641601563, 1559.37463378906, 2333.86547851563, 2112.19897460938, 5131.08935546875, 5671.3505859375, 3115.78466796875, 2336.2255859375, 2197.15649414063, 1073.25744628906, 1067.25158691406, 1029.73352050781, 1169.77075195313, 1487.70727539063, 1967.15612792969, 1458.63134765625, 2018.44934082031, 2272.38403320313, 1362.41247558594, 1066.52294921875, 1902.52453613281, 3767.75073242188, 4349.96630859375, 3587.845703125, 1167.42126464844, 1387.46203613281, 1277.7900390625, 907.506286621094, 1061.73120117188, 1034.74877929688, 4263.056640625 };

            var spectrumToDeconvolute = new MzmlMzSpectrum(mz1, intensities1, false);

            Console.WriteLine(string.Join(" , ", spectrumToDeconvolute.Deconvolute(4, new Tolerance(ToleranceUnit.Absolute, 0.01))));
        }

        [Test]
        public void SpectrumGetMassFirst()
        {
            Assert.AreEqual(328.73795, _mzSpectrumA.FirstX);
        }

        [Test]
        public void SpectrumGetMassRandom()
        {
            Assert.AreEqual(482.90393, _mzSpectrumA[6].Mz);
        }

        [Test]
        public void SpectrumContainsPeak()
        {
            Assert.IsTrue(_mzSpectrumA.Size > 0);
        }

        [Test]
        public void SpectrumContainsPeakInRange()
        {
            Assert.AreEqual(1, _mzSpectrumA.NumPeaksWithinRange(448.23987 - 0.001, 448.23987 + 0.001));
        }

        [Test]
        public void SpectrumContainsPeakInRangeEnd()
        {
            Assert.AreEqual(0, _mzSpectrumA.NumPeaksWithinRange(448.23987 - 0.001, 448.23987));
        }

        [Test]
        public void SpectrumContainsPeakInRangeStart()
        {
            Assert.AreEqual(1, _mzSpectrumA.NumPeaksWithinRange(448.23987, 448.23987 + 0.001));
        }

        [Test]
        public void SpectrumContainsPeakInRangeStartEnd()
        {
            Assert.AreEqual(0, _mzSpectrumA.NumPeaksWithinRange(448.23987, 448.23987));
        }

        [Test]
        public void SpectrumDoesntContainPeakInRange()
        {
            Assert.AreEqual(0, _mzSpectrumA.NumPeaksWithinRange(603.4243 - 0.001, 603.4243 + 0.001));
        }

        [Test]
        public void SpectrumMassRange()
        {
            MzRange range = new MzRange(328.73795, 723.35345);

            Assert.AreEqual(0, _mzSpectrumA.Range.Minimum - range.Minimum, 1e-9);
            Assert.AreEqual(0, _mzSpectrumA.Range.Maximum - range.Maximum, 1e-9);
        }

        [Test]
        public void SpectrumFilterCount()
        {
            var filteredMzSpectrum = _mzSpectrumA.FilterByY(28604417, 28604419);

            Assert.AreEqual(1, filteredMzSpectrum.Count());
        }

        [Test]
        public void SpectrumSelect()
        {
            MzSpectrum<MzmlPeak> v2 = _mzSpectrumA;
            ISpectrum<IPeak> v3 = v2;

            v3.Take(4);

            var v5 = v3.Select(b => b.X);
            Assert.AreEqual(328.73795, v5.First());

            var bn = v2[0];

            var bsrg = _mzSpectrumA[0];
        }

        [Test]
        public void FilterByNumberOfMostIntenseTest()
        {
            Assert.AreEqual(5, _mzSpectrumA.FilterByNumberOfMostIntense(5).Count());
        }

        [Test]
        public void GetBasePeak()
        {
            Assert.AreEqual(122781408.0, _mzSpectrumA.PeakWithHighestY.Intensity);
        }

        [Test]
        public void GetClosestPeak()
        {
            Assert.AreEqual(448.23987, _mzSpectrumA.GetClosestPeak(448).Mz);
            Assert.AreEqual(447.73849, _mzSpectrumA.GetClosestPeak(447.9).Mz);
        }

        [Test]
        public void Extract()
        {
            Assert.AreEqual(3, _mzSpectrumA.Extract(500, 600).Count());
        }

        [Test]
        public void CorrectOrder()
        {
            _mzSpectrumA = new MzmlMzSpectrum(new double[] { 5, 6, 7 }, new double[] { 1, 2, 3 }, false);
            Assert.IsTrue(_mzSpectrumA.FilterByNumberOfMostIntense(2).First().Mz < _mzSpectrumA.FilterByNumberOfMostIntense(2).ToList()[1].Mz);
        }

        [Test]
        public void TestFunctionToX()
        {
            _mzSpectrumA.ReplaceXbyApplyingFunction(b => -1);
            Assert.AreEqual(-1, _mzSpectrumA[0].X);
        }

        [Test]
        public void TestGetClosestPeakXValue()
        {
            Assert.AreEqual(447.73849, _mzSpectrumA.GetClosestPeakXvalue(447.73849));
            Assert.AreEqual(447.73849, _mzSpectrumA.GetClosestPeakXvalue(447));
            Assert.Throws<IndexOutOfRangeException>(() => { new MzmlMzSpectrum(new double[0], new double[0], false).GetClosestPeakXvalue(1); }, "No peaks in spectrum!");
        }

        //[Test]
        //public void IMzSpectrumTpeakTest()
        //{
        //    //double[] mz = { 328.73795, 329.23935, 447.73849, 448.23987, 482.23792, 482.57089, 482.90393, 500.95358, 501.28732, 501.62131, 611.99377, 612.32806, 612.66187, 722.85217, 723.35345 };
        //    //double[] intensities = { 81007096.0, 28604418.0, 78353512.0, 39291696.0, 122781408.0, 94147520.0, 44238040.0, 71198680.0, 54184096.0, 21975364.0, 44514172.0, 43061628.0, 23599424.0, 56022696.0, 41019144.0 };

        //    IMzSpectrum<MzmlPeak> ok = _mzSpectrumA;

        //    Assert.Greater(100, ok.ApplyFunctionToX(b => b / 10).LastX);

        //    Assert.AreEqual(1, ok.Extract(new DoubleRange(723, 1000)).Size);

        //    Assert.AreEqual(15, ok.Extract(double.MinValue, double.MaxValue).Size);

        //    Assert.AreEqual(0, ok.FilterByNumberOfMostIntense(0).Size);
        //    Assert.AreEqual(5, ok.FilterByNumberOfMostIntense(5).Size);
        //    Assert.AreEqual(15, ok.FilterByNumberOfMostIntense(15).Size);

        //    Assert.AreEqual(15, ok.FilterByY(new DoubleRange(double.MinValue, double.MaxValue)).Size);

        //    Assert.AreEqual(1, ok.FilterByY(39291695, 39291697).Size);

        //    Assert.AreEqual(2, ok.WithRangeRemoved(new DoubleRange(329, 723)).Size);

        //    Assert.AreEqual(15, ok.WithRangeRemoved(0, 1).Size);

        //    Assert.AreEqual(2, ok.WithRangesRemoved(new List<DoubleRange> { new DoubleRange(329, 400), new DoubleRange(400, 723) }).Size);
        //}

        [Test]
        public void TestNumPeaksWithinRange()
        {
            double[] xArray = { 1, 2, 3, 4, 5, 6, 7 };
            double[] yArray = { 1, 2, 1, 5, 1, 2, 1 };

            var thisSpectrum = new MzmlMzSpectrum(xArray, yArray, false);

            Assert.AreEqual(7, thisSpectrum.NumPeaksWithinRange(double.MinValue, double.MaxValue));

            Assert.AreEqual(6, thisSpectrum.NumPeaksWithinRange(1, 7));

            Assert.AreEqual(0, thisSpectrum.NumPeaksWithinRange(1, 1));

            Assert.AreEqual(1, thisSpectrum.NumPeaksWithinRange(1, 2));

            Assert.AreEqual(2, thisSpectrum.NumPeaksWithinRange(0.001, 2.999));

            Assert.AreEqual(1, thisSpectrum.NumPeaksWithinRange(0, 1.5));

            Assert.AreEqual(1, thisSpectrum.NumPeaksWithinRange(6.5, 8));

            Assert.AreEqual(2, thisSpectrum.NumPeaksWithinRange(3, 5));

            Assert.AreEqual(2, thisSpectrum.NumPeaksWithinRange(3.5, 5.5));

            Assert.AreEqual(1, thisSpectrum.NumPeaksWithinRange(7, 8));

            Assert.AreEqual(0, thisSpectrum.NumPeaksWithinRange(8, 9));

            Assert.AreEqual(0, thisSpectrum.NumPeaksWithinRange(-2, -1));

            Assert.AreEqual("[1 to 7] m/z (Peaks 7)", thisSpectrum.ToString());

            //Assert.AreEqual(7, thisSpectrum.FilterByNumberOfMostIntense(7).Size);
            //Assert.AreEqual(1, thisSpectrum.FilterByNumberOfMostIntense(1).Size);
            //Assert.AreEqual(4, thisSpectrum.FilterByNumberOfMostIntense(1).FirstX);

            //Assert.AreEqual(2, thisSpectrum.FilterByNumberOfMostIntense(3).FirstX);

            //Assert.AreEqual(0, thisSpectrum.FilterByNumberOfMostIntense(0).Size);

            //Assert.AreEqual(2, thisSpectrum.WithRangeRemoved(2, 6).Size);
            //Assert.AreEqual(0, thisSpectrum.WithRangeRemoved(0, 100).Size);

            //Assert.AreEqual(6, thisSpectrum.WithRangeRemoved(7, 100).Size);

            //Assert.AreEqual(1, thisSpectrum.WithRangeRemoved(new DoubleRange(double.MinValue, 6)).Size);

            List<DoubleRange> xRanges = new List<DoubleRange>
            {
                new DoubleRange(2, 5),
                new DoubleRange(3, 6)
            };
            //Assert.AreEqual(2, thisSpectrum.WithRangesRemoved(xRanges).Size);

            //Assert.AreEqual(3, thisSpectrum.Extract(new DoubleRange(4.5, 10)).Size);

            //Assert.AreEqual(2, thisSpectrum.FilterByY(new DoubleRange(1.5, 2.5)).Size);

            //Assert.AreEqual(3, thisSpectrum.FilterByY(1.5, double.MaxValue).Size);

            //Assert.AreEqual(2, thisSpectrum.ApplyFunctionToX(b => b * 2).FirstX);

            Assert.AreEqual(1, thisSpectrum.GetClosestPeak(-100).X);

            Assert.AreEqual(7, thisSpectrum.GetClosestPeak(6.6).X);

            Assert.AreEqual(7, thisSpectrum.GetClosestPeak(7).X);

            Assert.AreEqual(7, thisSpectrum.GetClosestPeak(8).X);

            IEnumerable hnm = thisSpectrum;

            double dudu = 0;
            foreach (var ikik in hnm)
            {
                dudu += ((Peak)ikik).X;
            }
            Assert.AreEqual(1 + 2 + 3 + 4 + 5 + 6 + 7, dudu);
        }

        #endregion Public Methods

    }
}