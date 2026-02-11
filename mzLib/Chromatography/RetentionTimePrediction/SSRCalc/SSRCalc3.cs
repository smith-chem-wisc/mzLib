using System.Text;

namespace Chromatography.RetentionTimePrediction.SSRCalc;

/**
/*
/* reference, O. V. Krokhin, R. Craig, V. Spicer, W. Ens, K. G. Standing, R. C. Beavis, J. A. Wilkins
/* An improved model for prediction of retention times of tryptic peptides in ion-pair reverse-phase HPLC:
/* its application to protein peptide mapping by off-line HPLC-MALDI MS
/* Molecular and Cellular Proteomics 2004 Sep;3(9):908-19.
/* URL, http://hs2.proteome.ca/SSRCalc/SSRCalc.html
/*
/*
/* These subroutines are based on web version SSRCalculator of the Copyright holder listed as in the following:
/*
/* Version 3.0   2005.02.28
/* Copyright (c) 2005 John Wilkins
/* Sequence Specific Retention Calculator
/* Authors: Oleg Krokhin, Vic Spicer, John Cortens
 */

/* Translated from perl to C, Ted Holzman FHCRC, 6/2006  */
/* Retranslated from C to Java, Ted Holzman FHCRC 7/2006 */
/* Translated from Java to C#, Brendan MacLean UW 10/2008 */
/* NB: This is a version 0.1 direct translation.
/*     An attempt has been made to keep function names, variable names, and algorithms
/*     as close as possible to the original perl.
 */

public class SSRCalc3
{
    public const String VERSION = "Krokhin,3.0";

    private static readonly CLUSTCOMB_List CLUSTCOMB = new CLUSTCOMB_List();
    private static readonly Dictionary<string, double> HlxScore4 = new Dictionary<string, double>();
    private static readonly Dictionary<string, double> HlxScore5 = new Dictionary<string, double>();
    private static readonly Dictionary<string, double> HlxScore6 = new Dictionary<string, double>();
    private static readonly int[] EMap = new int[128];

    private sealed class CLUSTCOMB_List : List<KeyValuePair<string, double>>
    {
        public void Add(string pattern, double value)
        {
            Add(new KeyValuePair<string, double>(pattern, value));
        }
    }

    static SSRCalc3()
    {
        // ReSharper disable NonLocalizedString
        CLUSTCOMB.Add("0110", 0.3);
        CLUSTCOMB.Add("0150", 0.4);
        CLUSTCOMB.Add("0510", 0.4);
        CLUSTCOMB.Add("0550", 1.3);
        CLUSTCOMB.Add("01110", 0.5);
        CLUSTCOMB.Add("01150", 0.7);
        CLUSTCOMB.Add("01510", 0.7);
        CLUSTCOMB.Add("01550", 2.1);
        CLUSTCOMB.Add("05110", 0.7);
        CLUSTCOMB.Add("05150", 2.1);
        CLUSTCOMB.Add("05510", 2.1);
        CLUSTCOMB.Add("05550", 2.8);
        CLUSTCOMB.Add("011110", 0.7);
        CLUSTCOMB.Add("011150", 0.9);
        CLUSTCOMB.Add("011510", 0.9);
        CLUSTCOMB.Add("011550", 2.2);
        CLUSTCOMB.Add("015110", 0.9);
        CLUSTCOMB.Add("015150", 2.2);
        CLUSTCOMB.Add("015510", 0.9);
        CLUSTCOMB.Add("015550", 3.0);
        CLUSTCOMB.Add("051110", 0.9);
        CLUSTCOMB.Add("051150", 2.2);
        CLUSTCOMB.Add("051510", 2.2);
        CLUSTCOMB.Add("051550", 3.0);
        CLUSTCOMB.Add("055110", 2.2);
        CLUSTCOMB.Add("055150", 3.0);
        CLUSTCOMB.Add("055510", 3.0);
        CLUSTCOMB.Add("055550", 3.5);
        CLUSTCOMB.Add("0111110", 0.9);
        CLUSTCOMB.Add("0111150", 1.0);
        CLUSTCOMB.Add("0111510", 1.0);
        CLUSTCOMB.Add("0111550", 2.3);
        CLUSTCOMB.Add("0115110", 1.0);
        CLUSTCOMB.Add("0115150", 2.3);
        CLUSTCOMB.Add("0115510", 2.3);
        CLUSTCOMB.Add("0115550", 3.1);
        CLUSTCOMB.Add("0151110", 1.0);
        CLUSTCOMB.Add("0151150", 2.3);
        CLUSTCOMB.Add("0151510", 2.3);
        CLUSTCOMB.Add("0151550", 3.1);
        CLUSTCOMB.Add("0155110", 2.3);
        CLUSTCOMB.Add("0155150", 3.1);
        CLUSTCOMB.Add("0155510", 3.1);
        CLUSTCOMB.Add("0155550", 3.6);
        CLUSTCOMB.Add("0511110", 1.0);
        CLUSTCOMB.Add("0511150", 2.3);
        CLUSTCOMB.Add("0511510", 2.3);
        CLUSTCOMB.Add("0511550", 3.1);
        CLUSTCOMB.Add("0515110", 3.6);
        CLUSTCOMB.Add("0515150", 2.3);
        CLUSTCOMB.Add("0515510", 3.1);
        CLUSTCOMB.Add("0515550", 3.6);
        CLUSTCOMB.Add("0551110", 2.3);
        CLUSTCOMB.Add("0551150", 3.1);
        CLUSTCOMB.Add("0551510", 3.1);
        CLUSTCOMB.Add("0551550", 3.6);
        CLUSTCOMB.Add("0555110", 3.1);
        CLUSTCOMB.Add("0555150", 3.6);
        CLUSTCOMB.Add("0555510", 3.6);
        CLUSTCOMB.Add("0555550", 4.0);
        CLUSTCOMB.Add("01111110", 1.1);
        CLUSTCOMB.Add("01111150", 1.7);
        CLUSTCOMB.Add("01111510", 1.7);
        CLUSTCOMB.Add("01111550", 2.5);
        CLUSTCOMB.Add("01115110", 1.7);
        CLUSTCOMB.Add("01115150", 2.5);
        CLUSTCOMB.Add("01115510", 2.5);
        CLUSTCOMB.Add("01115550", 3.3);
        CLUSTCOMB.Add("01151110", 1.7);
        CLUSTCOMB.Add("01151150", 2.5);
        CLUSTCOMB.Add("01151510", 2.5);
        CLUSTCOMB.Add("01151550", 3.3);
        CLUSTCOMB.Add("01155110", 2.5);
        CLUSTCOMB.Add("01155150", 3.3);
        CLUSTCOMB.Add("01155510", 3.3);
        CLUSTCOMB.Add("01155550", 3.7);
        CLUSTCOMB.Add("01511110", 1.7);
        CLUSTCOMB.Add("01511150", 2.5);
        CLUSTCOMB.Add("01511510", 2.5);
        CLUSTCOMB.Add("01511550", 3.3);
        CLUSTCOMB.Add("01515110", 2.5);
        CLUSTCOMB.Add("01515150", 3.3);
        CLUSTCOMB.Add("01515510", 3.3);
        CLUSTCOMB.Add("01515550", 3.7);
        CLUSTCOMB.Add("01551110", 2.5);
        CLUSTCOMB.Add("01551150", 3.3);
        CLUSTCOMB.Add("01551510", 3.3);
        CLUSTCOMB.Add("01551550", 3.7);
        CLUSTCOMB.Add("01555110", 3.3);
        CLUSTCOMB.Add("01555150", 3.7);
        CLUSTCOMB.Add("01555510", 3.7);
        CLUSTCOMB.Add("01555550", 4.1);
        CLUSTCOMB.Add("05111110", 1.7);
        CLUSTCOMB.Add("05111150", 2.5);
        CLUSTCOMB.Add("05111510", 2.5);
        CLUSTCOMB.Add("05111550", 3.3);
        CLUSTCOMB.Add("05115110", 2.5);
        CLUSTCOMB.Add("05115150", 3.3);
        CLUSTCOMB.Add("05115510", 3.3);
        CLUSTCOMB.Add("05115550", 3.7);
        CLUSTCOMB.Add("05151110", 2.5);
        CLUSTCOMB.Add("05151150", 3.3);
        CLUSTCOMB.Add("05151510", 3.3);
        CLUSTCOMB.Add("05151550", 3.7);
        CLUSTCOMB.Add("05155110", 3.3);
        CLUSTCOMB.Add("05155150", 3.7);
        CLUSTCOMB.Add("05155510", 3.7);
        CLUSTCOMB.Add("05155550", 4.1);
        CLUSTCOMB.Add("05511110", 2.5);
        CLUSTCOMB.Add("05511150", 3.3);
        CLUSTCOMB.Add("05511510", 3.3);
        CLUSTCOMB.Add("05511550", 3.7);
        CLUSTCOMB.Add("05515110", 3.3);
        CLUSTCOMB.Add("05515150", 3.7);
        CLUSTCOMB.Add("05515510", 3.7);
        CLUSTCOMB.Add("05515550", 4.1);
        CLUSTCOMB.Add("05551110", 3.3);
        CLUSTCOMB.Add("05551150", 3.7);
        CLUSTCOMB.Add("05551510", 3.7);
        CLUSTCOMB.Add("05551550", 4.1);
        CLUSTCOMB.Add("05555110", 3.7);
        CLUSTCOMB.Add("05555150", 4.1);
        CLUSTCOMB.Add("05555510", 4.1);
        CLUSTCOMB.Add("05555550", 4.5);

        HlxScore4.Add("XXUX", 0.8);
        HlxScore4.Add("XZOX", 0.8);
        HlxScore4.Add("XUXX", 0.8);
        HlxScore4.Add("XXOX", 0.7);
        HlxScore4.Add("XOXX", 0.7);
        HlxScore4.Add("XZUX", 0.7);
        HlxScore4.Add("XXOZ", 0.7);
        HlxScore4.Add("ZXOX", 0.7);
        HlxScore4.Add("XOZZ", 0.7);
        HlxScore4.Add("ZOXX", 0.7);
        HlxScore4.Add("ZOZX", 0.7);
        HlxScore4.Add("ZUXX", 0.7);
        HlxScore4.Add("ZXUX", 0.5);
        HlxScore4.Add("XOZX", 0.5);
        HlxScore4.Add("XZOZ", 0.5);
        HlxScore4.Add("XUZX", 0.5);
        HlxScore4.Add("ZZOX", 0.2);
        HlxScore4.Add("ZXOZ", 0.2);
        HlxScore4.Add("ZOXZ", 0.2);
        HlxScore4.Add("XOXZ", 0.2);
        HlxScore4.Add("ZZUZ", 0.2);
        HlxScore4.Add("XUXZ", 0.2);
        HlxScore4.Add("ZUXZ", 0.2);
        HlxScore4.Add("XZUZ", 0.2);
        HlxScore4.Add("XUZZ", 0.2);
        HlxScore4.Add("ZXUZ", 0.2);
        HlxScore4.Add("ZOZZ", 0.2);
        HlxScore4.Add("ZZOZ", 0.2);
        HlxScore4.Add("ZZUX", 0.2);
        HlxScore4.Add("ZUZX", 0.2);
        HlxScore4.Add("XXUZ", 0.2);
        HlxScore4.Add("ZUZZ", 0.2);

        HlxScore5.Add("XXOXX", 3.75);
        HlxScore5.Add("XXOXZ", 3.75);
        HlxScore5.Add("XXOZX", 3.75);
        HlxScore5.Add("XZOXX", 3.75);
        HlxScore5.Add("ZXOXX", 3.75);
        HlxScore5.Add("XXOZZ", 2.7);
        HlxScore5.Add("XZOXZ", 2.7);
        HlxScore5.Add("XZOZX", 2.7);
        HlxScore5.Add("ZXOXZ", 2.7);
        HlxScore5.Add("ZXOZX", 2.7);
        HlxScore5.Add("ZZOXX", 2.7);
        HlxScore5.Add("ZXOZZ", 1.3);
        HlxScore5.Add("XZOZZ", 1.3);
        HlxScore5.Add("ZZOXZ", 1.3);
        HlxScore5.Add("ZZOZX", 1.3);
        HlxScore5.Add("ZZOZZ", 1.3);
        HlxScore5.Add("XXUXX", 3.75);
        HlxScore5.Add("XXUXZ", 3.75);
        HlxScore5.Add("XXUZX", 3.75);
        HlxScore5.Add("XZUXX", 3.75);
        HlxScore5.Add("ZXUXX", 3.75);
        HlxScore5.Add("XXUZZ", 1.1);
        HlxScore5.Add("XZUXZ", 1.1);
        HlxScore5.Add("XZUZX", 1.1);
        HlxScore5.Add("ZXUZX", 1.1);
        HlxScore5.Add("ZXUXZ", 1.1);
        HlxScore5.Add("ZZUXX", 1.1);
        HlxScore5.Add("XZUZZ", 1.3);
        HlxScore5.Add("ZXUZZ", 1.3);
        HlxScore5.Add("ZZUXZ", 1.3);
        HlxScore5.Add("ZZUZX", 1.3);
        HlxScore5.Add("ZZUZZ", 1.3);
        HlxScore5.Add("XXOOX", 1.25);
        HlxScore5.Add("ZXOOX", 1.25);
        HlxScore5.Add("XZOOX", 1.25);
        HlxScore5.Add("XOOXX", 1.25);
        HlxScore5.Add("XOOXZ", 1.25);
        HlxScore5.Add("XOOZX", 1.25);
        HlxScore5.Add("XXOOZ", 1.25);
        HlxScore5.Add("ZXOOZ", 1.25);
        HlxScore5.Add("XZOOZ", 1.25);
        HlxScore5.Add("ZZOOX", 1.25);
        HlxScore5.Add("ZZOOZ", 1.25);
        HlxScore5.Add("ZOOXX", 1.25);
        HlxScore5.Add("ZOOXZ", 1.25);
        HlxScore5.Add("ZOOZX", 1.25);
        HlxScore5.Add("XOOZZ", 1.25);
        HlxScore5.Add("ZOOZZ", 1.25);
        HlxScore5.Add("XXOUX", 1.25);
        HlxScore5.Add("ZXOUX", 1.25);
        HlxScore5.Add("XXUOX", 1.25);
        HlxScore5.Add("ZXUOX", 1.25);
        HlxScore5.Add("XOUXX", 1.25);
        HlxScore5.Add("XOUXZ", 1.25);
        HlxScore5.Add("XUOXX", 1.25);
        HlxScore5.Add("XUOXZ", 1.25);
        HlxScore5.Add("XXOUZ", 0.75);
        HlxScore5.Add("ZXOUZ", 0.75);
        HlxScore5.Add("XZOUX", 0.75);
        HlxScore5.Add("XZOUZ", 0.75);
        HlxScore5.Add("ZZOUX", 0.75);
        HlxScore5.Add("ZZOUZ", 0.75);
        HlxScore5.Add("XXUOZ", 0.75);
        HlxScore5.Add("ZXUOZ", 0.75);
        HlxScore5.Add("XZUOX", 0.75);
        HlxScore5.Add("XZUOZ", 0.75);
        HlxScore5.Add("ZZUOX", 0.75);
        HlxScore5.Add("ZZUOZ", 0.75);
        HlxScore5.Add("ZOUXX", 0.75);
        HlxScore5.Add("ZOUXZ", 0.75);
        HlxScore5.Add("XOUZX", 0.75);
        HlxScore5.Add("ZOUZX", 0.75);
        HlxScore5.Add("XOUZZ", 0.75);
        HlxScore5.Add("ZOUZZ", 0.75);
        HlxScore5.Add("ZUOXX", 0.75);
        HlxScore5.Add("ZUOXZ", 0.75);
        HlxScore5.Add("XUOZX", 0.75);
        HlxScore5.Add("ZUOZX", 0.75);
        HlxScore5.Add("XUOZZ", 0.75);
        HlxScore5.Add("ZUOZZ", 0.75);
        HlxScore5.Add("XUUXX", 1.25);
        HlxScore5.Add("XXUUX", 1.25);
        HlxScore5.Add("XXUUZ", 0.6);
        HlxScore5.Add("ZXUUX", 0.6);
        HlxScore5.Add("ZXUUZ", 0.6);
        HlxScore5.Add("XZUUX", 0.6);
        HlxScore5.Add("XZUUZ", 0.6);
        HlxScore5.Add("ZZUUX", 0.6);
        HlxScore5.Add("ZZUUZ", 0.6);
        HlxScore5.Add("ZUUXX", 0.6);
        HlxScore5.Add("XUUXZ", 0.6);
        HlxScore5.Add("ZUUXZ", 0.6);
        HlxScore5.Add("XUUZX", 0.6);
        HlxScore5.Add("ZUUZX", 0.6);
        HlxScore5.Add("XUUZZ", 0.6);
        HlxScore5.Add("ZUUZZ", 0.6);

        HlxScore6.Add("XXOOXX", 3.0);
        HlxScore6.Add("XXOOXZ", 3.0);
        HlxScore6.Add("ZXOOXX", 3.0);
        HlxScore6.Add("ZXOOXZ", 3.0);
        HlxScore6.Add("XXOUXX", 3.0);
        HlxScore6.Add("XXOUXZ", 3.0);
        HlxScore6.Add("XXUOXX", 3.0);
        HlxScore6.Add("XXUOXZ", 3.0);
        HlxScore6.Add("ZXUOXX", 3.0);
        HlxScore6.Add("ZXOUXX", 3.0);
        HlxScore6.Add("XXOOZX", 1.6);
        HlxScore6.Add("XXOOZZ", 1.6);
        HlxScore6.Add("XZOOXX", 1.6);
        HlxScore6.Add("XZOOXZ", 1.6);
        HlxScore6.Add("XZOOZX", 1.6);
        HlxScore6.Add("XZOOZZ", 1.6);
        HlxScore6.Add("ZXOOZX", 1.6);
        HlxScore6.Add("ZXOOZZ", 1.6);
        HlxScore6.Add("ZZOOXX", 1.6);
        HlxScore6.Add("ZZOOXZ", 1.6);
        HlxScore6.Add("ZXOUXZ", 1.6);
        HlxScore6.Add("XZUOXX", 1.6);
        HlxScore6.Add("ZXUOXZ", 1.6);
        HlxScore6.Add("ZZOOZX", 1.5);
        HlxScore6.Add("ZZOOZZ", 1.5);
        HlxScore6.Add("XXOUZX", 1.5);
        HlxScore6.Add("XXOUZZ", 1.5);
        HlxScore6.Add("XZOUXX", 1.5);
        HlxScore6.Add("XZOUXZ", 1.5);
        HlxScore6.Add("ZXOUZX", 1.5);
        HlxScore6.Add("ZXOUZZ", 1.5);
        HlxScore6.Add("ZZOUXX", 1.5);
        HlxScore6.Add("ZZOUXZ", 1.5);
        HlxScore6.Add("XXUOZX", 1.5);
        HlxScore6.Add("XXUOZZ", 1.5);
        HlxScore6.Add("XZUOXZ", 1.5);
        HlxScore6.Add("ZXUOZX", 1.5);
        HlxScore6.Add("ZXUOZZ", 1.5);
        HlxScore6.Add("ZZUOXX", 1.5);
        HlxScore6.Add("ZZUOXZ", 1.5);
        HlxScore6.Add("ZZUOZX", 1.25);
        HlxScore6.Add("ZZUOZZ", 1.25);
        HlxScore6.Add("ZZOUZX", 1.25);
        HlxScore6.Add("ZZOUZZ", 1.25);
        HlxScore6.Add("XZOUZX", 1.25);
        HlxScore6.Add("XZOUZZ", 1.25);
        HlxScore6.Add("XZUOZX", 1.25);
        HlxScore6.Add("XZUOZZ", 1.25);
        HlxScore6.Add("XXUUXX", 1.25);
        HlxScore6.Add("XXUUXZ", 1.25);
        HlxScore6.Add("ZXUUXX", 1.25);
        HlxScore6.Add("XXUUZX", 1.25);
        HlxScore6.Add("XXUUZZ", 1.25);
        HlxScore6.Add("XZUUXX", 1.25);
        HlxScore6.Add("XZUUXZ", 1.25);
        HlxScore6.Add("XZUUZX", 0.75);
        HlxScore6.Add("XZUUZZ", 0.75);
        HlxScore6.Add("ZXUUXZ", 1.25);
        HlxScore6.Add("ZXUUZX", 1.25);
        HlxScore6.Add("ZXUUZZ", 1.25);
        HlxScore6.Add("ZZUUXX", 1.25);
        HlxScore6.Add("ZZUUXZ", 1.25);
        HlxScore6.Add("ZZUUZX", 0.75);
        HlxScore6.Add("ZZUUZZ", 0.75);
        // ReSharper restore NonLocalizedString

        for (int i = 0; i < EMap.Length; i++)
        {
            EMap[i] = -1;
        }

        EMap['K'] = 0;
        EMap['R'] = 1;
        EMap['H'] = 2;
        EMap['D'] = 3;
        EMap['E'] = 4;
        EMap['C'] = 5;
        EMap['Y'] = 6;
    }

    public enum Column { A300, A100 }
    public AAParams[] AAPARAMS = new AAParams[128];

    public SSRCalc3(string name, Column column)
    {
        Name = name;
        AAParams NULLPARAM = new AAParams(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

        for (int i = 0; i < AAPARAMS.Length; i++)
        {
            AAPARAMS[i] = NULLPARAM;
        }

        switch (column)
        {
            case Column.A300:
                A300Column();
                break;
            case Column.A100:
                A100Column();
                break;
        }
    }

    public string Name { get; private set; }

    private void A300Column()
    {
        AAPARAMS['A'] = new AAParams(01.10, 00.35, 00.50, 00.80, -0.10, 00.80, -0.30, 00.10, 00.80, -0.50, 00.00, 071.0370, 3.55, 7.59, 00.00, 1.0, 1.2);
        AAPARAMS['C'] = new AAParams(00.45, 00.90, 00.20, -0.80, -0.50, 00.50, 00.40, 00.00, -0.80, -0.50, 00.00, 103.0090, 3.55, 7.50, 00.00, 0.0, 1.0);
        AAPARAMS['D'] = new AAParams(00.15, 00.50, 00.40, -0.50, -0.50, 00.30, 00.30, 00.70, -0.50, -0.50, 00.00, 115.0270, 4.55, 7.50, 04.05, 0.0, 1.1);
        AAPARAMS['E'] = new AAParams(00.95, 01.00, 00.00, 00.00, -0.10, 00.50, 00.10, 00.00, 00.00, -0.10, 00.00, 129.0430, 4.75, 7.70, 04.45, 0.0, 1.1);
        AAPARAMS['F'] = new AAParams(10.90, 07.50, 09.50, 10.50, 10.30, 11.10, 08.10, 09.50, 10.50, 10.30, -0.10, 147.0638, 3.55, 7.50, 00.00, 0.5, 1.0);
        AAPARAMS['G'] = new AAParams(-0.35, 00.20, 00.15, -0.90, -0.70, 00.00, 00.00, 00.10, -0.90, -0.70, 00.00, 057.0210, 3.55, 7.50, 00.00, 0.0, 0.3);
        AAPARAMS['H'] = new AAParams(-1.45, -0.10, -0.20, -1.30, -1.70, -1.00, 00.10, -0.20, -1.30, -1.70, 00.00, 137.0590, 3.55, 7.50, 05.98, 0.0, 0.6);
        AAPARAMS['I'] = new AAParams(08.00, 05.20, 06.60, 08.40, 07.70, 07.70, 05.00, 06.80, 08.40, 07.70, 00.15, 113.0840, 3.55, 7.50, 00.00, 3.5, 1.4);
        AAPARAMS['K'] = new AAParams(-2.05, -0.60, -1.50, -1.90, -1.45, -0.20, -1.40, -1.30, -2.20, -1.45, 00.00, 128.0950, 3.55, 7.50, 10.00, 0.0, 1.0);
        AAPARAMS['L'] = new AAParams(09.30, 05.55, 07.40, 09.60, 09.30, 09.20, 06.00, 07.90, 09.60, 08.70, 00.30, 113.0840, 3.55, 7.50, 00.00, 1.6, 1.6);
        AAPARAMS['M'] = new AAParams(06.20, 04.40, 05.70, 05.80, 06.00, 06.20, 05.00, 05.70, 05.80, 06.00, 00.00, 131.0400, 3.55, 7.00, 00.00, 1.8, 1.0);
        AAPARAMS['N'] = new AAParams(-0.85, 00.20, -0.20, -1.20, -1.10, -0.85, 00.20, -0.20, -1.20, -1.10, 00.00, 114.0430, 3.55, 7.50, 00.00, 0.0, 0.4);
        AAPARAMS['P'] = new AAParams(02.10, 02.10, 02.10, 00.20, 02.10, 03.00, 01.00, 01.50, 00.20, 02.10, 00.00, 097.0530, 3.55, 8.36, 00.00, 0.0, 0.3);
        AAPARAMS['Q'] = new AAParams(-0.40, -0.70, -0.20, -0.90, -1.10, -0.40, -0.80, -0.20, -0.90, -1.10, 00.00, 128.0590, 3.55, 7.50, 00.00, 0.0, 1.0);
        AAPARAMS['R'] = new AAParams(-1.40, 00.50, -1.10, -1.30, -1.10, -0.20, 00.50, -1.10, -1.20, -1.10, 00.00, 156.1010, 3.55, 7.50, 12.00, 0.0, 1.0);
        AAPARAMS['S'] = new AAParams(-0.15, 00.80, -0.10, -0.80, -1.20, -0.50, 00.40, 00.10, -0.80, -1.20, 00.00, 087.0320, 3.55, 6.93, 00.00, 0.0, 1.0);
        AAPARAMS['T'] = new AAParams(00.65, 00.80, 00.60, 00.40, 00.00, 00.60, 00.80, 00.40, 00.40, 00.00, 00.00, 101.0480, 3.55, 6.82, 00.00, 0.0, 1.0);
        AAPARAMS['V'] = new AAParams(05.00, 02.90, 03.40, 05.00, 04.20, 05.10, 02.70, 03.40, 05.00, 04.20, -0.30, 099.0680, 3.55, 7.44, 00.00, 1.4, 1.2);
        AAPARAMS['W'] = new AAParams(12.25, 11.10, 11.80, 11.00, 12.10, 12.40, 11.60, 11.80, 11.00, 12.10, 00.15, 186.0790, 3.55, 7.50, 00.00, 1.6, 1.0);
        AAPARAMS['Y'] = new AAParams(04.85, 03.70, 04.50, 04.00, 04.40, 05.10, 04.20, 04.50, 04.00, 04.40, -0.20, 163.0630, 3.55, 7.50, 10.00, 0.2, 1.0);
        AAPARAMS['B'] = new AAParams(00.15, 00.50, 00.40, -0.50, -0.50, 00.30, 00.30, 00.70, -0.50, -0.50, 00.00, 115.0270, 4.55, 7.50, 04.05, 0.0, 1.1);
        AAPARAMS['X'] = new AAParams(00.00, 00.00, 00.00, 00.00, 00.00, 00.00, 00.00, 00.00, 00.00, 00.00, 00.00, 000.0000, 0.00, 0.00, 00.00, 0.0, 1.0);
        AAPARAMS['Z'] = new AAParams(00.95, 01.00, 00.00, 00.00, -0.10, 00.50, 00.10, 00.00, 00.00, -0.10, 00.00, 129.0430, 4.75, 7.70, 04.45, 0.0, 1.1);
    }

    private void A100Column()
    {
        AAPARAMS['A'] = new AAParams(01.02, -0.35, 00.35, 01.02, -0.20, 00.50, -0.05, 00.10, 00.50, -0.30, 00.00, 071.0370, 3.55, 7.59, 00.00, 1.0, 1.2);
        AAPARAMS['C'] = new AAParams(00.10, 00.40, 00.20, 00.10, -0.40, 00.60, 00.60, 01.00, 00.60, -0.50, 00.00, 103.0090, 3.55, 7.50, 00.00, 0.0, 1.0);
        AAPARAMS['D'] = new AAParams(00.15, 00.90, 00.60, 00.15, -0.40, 00.60, 00.30, 00.20, 00.60, -0.50, 00.00, 115.0270, 4.55, 7.50, 04.05, 0.0, 1.1);
        AAPARAMS['E'] = new AAParams(01.00, 01.00, -0.20, 01.00, -0.10, 00.70, 00.45, 00.50, 00.00, 00.25, 00.00, 129.0430, 4.75, 7.70, 04.45, 0.0, 1.1);
        AAPARAMS['F'] = new AAParams(11.67, 07.60, 09.70, 11.67, 11.50, 11.30, 08.40, 10.00, 11.30, 10.85, -0.10, 147.0638, 3.55, 7.50, 00.00, 0.5, 1.0);
        AAPARAMS['G'] = new AAParams(-0.35, 00.15, 00.15, -0.35, -0.40, 00.00, 00.15, 00.20, 00.00, -0.70, 00.00, 057.0210, 3.55, 7.50, 00.00, 0.0, 0.3);
        AAPARAMS['H'] = new AAParams(-3.00, -1.40, -1.00, -3.00, -1.90, -1.30, -1.30, -1.10, -1.30, -1.70, 00.00, 137.0590, 3.55, 7.50, 05.98, 0.0, 0.6);
        AAPARAMS['I'] = new AAParams(07.96, 04.95, 06.30, 07.96, 06.60, 07.25, 04.50, 06.50, 07.25, 07.20, 00.15, 113.0840, 3.55, 7.50, 00.00, 3.5, 1.4);
        AAPARAMS['K'] = new AAParams(-3.40, -1.85, -2.30, -2.10, -2.10, -1.75, -1.50, -1.75, -2.30, -2.50, 00.00, 128.0950, 3.55, 7.50, 10.00, 0.0, 1.0);
        AAPARAMS['L'] = new AAParams(09.40, 05.57, 07.40, 09.40, 09.30, 08.70, 05.50, 07.70, 08.70, 08.50, 00.30, 113.0840, 3.55, 7.50, 00.00, 1.6, 1.6);
        AAPARAMS['M'] = new AAParams(06.27, 05.20, 05.70, 06.27, 05.80, 06.25, 04.20, 05.70, 06.25, 05.60, 00.00, 131.0400, 3.55, 7.00, 00.00, 1.8, 1.0);
        AAPARAMS['N'] = new AAParams(-0.95, 01.20, -0.10, -0.95, -1.30, -0.65, 00.40, -0.05, -0.65, -1.20, 00.00, 114.0430, 3.55, 7.50, 00.00, 0.0, 0.4);
        AAPARAMS['P'] = new AAParams(01.85, 01.70, 01.75, 01.85, 01.20, 02.50, 01.70, 02.10, 02.50, 01.90, 00.00, 097.0530, 3.55, 8.36, 00.00, 0.0, 0.3);
        AAPARAMS['Q'] = new AAParams(-0.60, -0.50, -0.20, -0.60, -1.10, -0.40, -0.20, -0.70, -0.40, -1.30, 00.00, 128.0590, 3.55, 7.50, 00.00, 0.0, 1.0);
        AAPARAMS['R'] = new AAParams(-2.55, -1.40, -1.50, -1.10, -1.30, -1.00, 00.40, -1.00, -1.10, -1.90, 00.00, 156.1010, 3.55, 7.50, 12.00, 0.0, 1.0);
        AAPARAMS['S'] = new AAParams(-0.14, 01.10, -0.10, -0.14, -1.00, -0.40, 00.20, -0.30, -0.40, -1.20, 00.00, 087.0320, 3.55, 6.93, 00.00, 0.0, 1.0);
        AAPARAMS['T'] = new AAParams(00.64, 00.95, 00.60, 00.64, -0.10, 00.40, 00.30, 00.40, 00.40, -0.50, 00.00, 101.0480, 3.55, 6.82, 00.00, 0.0, 1.0);
        AAPARAMS['V'] = new AAParams(04.68, 02.10, 03.40, 04.68, 03.90, 04.40, 02.10, 03.00, 04.40, 04.40, -0.30, 099.0680, 3.55, 7.44, 00.00, 1.4, 1.2);
        AAPARAMS['W'] = new AAParams(13.35, 11.50, 11.80, 13.35, 13.00, 13.90, 11.80, 13.00, 13.90, 12.90, 00.15, 186.0790, 3.55, 7.50, 00.00, 1.6, 1.0);
        AAPARAMS['Y'] = new AAParams(05.35, 04.30, 05.10, 05.35, 05.00, 05.70, 05.00, 05.40, 05.70, 05.30, -0.20, 163.0630, 3.55, 7.50, 10.00, 0.2, 1.0);
        AAPARAMS['B'] = new AAParams(00.15, 00.50, 00.40, -0.50, -0.50, 00.30, 00.30, 00.70, -0.50, -0.50, 00.00, 115.0270, 4.55, 7.50, 04.05, 0.0, 1.1);
        AAPARAMS['X'] = new AAParams(00.00, 00.00, 00.00, 00.00, 00.00, 00.00, 00.00, 00.00, 00.00, 00.00, 00.00, 000.0000, 0.00, 0.00, 00.00, 0.0, 1.0);
        AAPARAMS['Z'] = new AAParams(00.95, 01.00, 00.00, 00.00, -0.10, 00.50, 00.10, 00.00, 00.00, -0.10, 00.00, 129.0430, 4.75, 7.70, 04.45, 0.0, 1.1);
    }

    public int NOELECTRIC { get; set; }
    public int NOCLUSTER { get; set; }
    public int NODIGEST { get; set; }
    public int NOSMALL { get; set; }
    public int NOHELIX1 { get; set; }
    public int NOHELIX2 { get; set; }
    public int NOEHEL { get; set; }

    private const bool DUPLICATE_ORIGINAL_CODE = true;
    //Speaking with the developers, it
    // was determined that it ought not to have been commented out.  So --
    // ALGORITHM_VERSION may be used to choose the older or newer code
    private const int ALGORITHM_VERSION = 3;
    private const int LPLim = 20;
    private const int SPLim = 8;
    private const double LPSFac = 0.0270;
    private const double SPSFac = -0.055;
    private const double UDF21 = 0.0, UDF22 = 0.0;
    private const double UDF31 = 1.0, UDF32 = 0.0;
    private const double SUMSCALE1 = 0.27, SUMSCALE2 = 0.33, SUMSCALE3 = 0.38, SUMSCALE4 = 0.447;
    private const double KSCALE = 0.4;
    private const double Z01 = -0.03, Z02 = 0.60, NDELTAWT = 0.8;
    private const double Z03 = 0.00, Z04 = 0.00, PDELTAWT = 1.0;
    private const double PPSCORE = 1.2, PPPSCORE = 3.5, PPPPSCORE = 5.0;
    private const double HELIX1SCALE = 1.6, HELIX2SCALE = 0.255;
    private const int HISC = 0;
    private const int GSC = 1;

    public double UnknownScore => 0;

    public double ScoreSequence(string baseSequence)
    {
        var seq = baseSequence;
        double tsum3 = 0.0;
        int sze = seq.Length;
        if (sze < 4) return tsum3;

        if (sze < 10)
        {
            tsum3 = AAPARAMS[seq[0]].RC1S + AAPARAMS[seq[1]].RC2S + AAPARAMS[seq[sze - 1]].RCNS + AAPARAMS[seq[sze - 2]].RCN2S;
            for (int i = 2; i < sze - 2; i++) tsum3 += AAPARAMS[seq[i]].RCS;
        }
        else
        {
            tsum3 = AAPARAMS[seq[0]].RC1 + AAPARAMS[seq[1]].RC2 + AAPARAMS[seq[sze - 1]].RCN + AAPARAMS[seq[sze - 2]].RCN2;
            for (int i = 2; i < sze - 2; i++) tsum3 += AAPARAMS[seq[i]].RC;
        }

        tsum3 += Smallness(sze, tsum3);
        tsum3 -= Undigested(seq);
        tsum3 -= Clusterness(seq);
        tsum3 -= Proline(seq);
        tsum3 *= Length_scale(sze);
        if (tsum3 >= 20 && tsum3 < 30) tsum3 -= ((tsum3 - 18) * SUMSCALE1);
        if (tsum3 >= 30 && tsum3 < 40) tsum3 -= ((tsum3 - 18) * SUMSCALE2);
        if (tsum3 >= 40 && tsum3 < 50) tsum3 -= ((tsum3 - 18) * SUMSCALE3);
        if (tsum3 >= 50) tsum3 -= ((tsum3 - 18) * SUMSCALE4);
        tsum3 += NewIso(seq, tsum3);
        tsum3 += Helicity1(seq);
        tsum3 += Helicity2(seq);
        tsum3 += Helectric(seq);
        return tsum3;
    }

    private double Smallness(int sqlen, double tsum)
    {
        if (NOSMALL == 1) return 0.0;
        if (sqlen < 20 && (tsum / sqlen) < 0.9) return 3.5 * (0.9 - (tsum / sqlen));
        if (sqlen < 15 && (tsum / sqlen) > 2.8) return 2.6 * ((tsum / sqlen) - 2.8);
        return 0.0;
    }

    private double Undigested(string sq)
    {
        if (NODIGEST == 1) return 0.0;
        char op1, op2;
        int xx = sq.Length - 1;
        char re = sq[xx];
        double csum = 0.0;
        if (re == 'R' || re == 'K' || re == 'H')
        {
            op1 = sq[xx - 1]; op2 = sq[xx - 2];
            csum = UDF21 * AAPARAMS[op1].UndKRH + UDF22 * AAPARAMS[op2].UndKRH;
        }
        for (int dd = 0; dd < xx; dd++)
        {
            re = sq[dd];
            if (re == 'K' || re == 'R' || re == 'H')
            {
                char op3, op4;
                op1 = op2 = op3 = op4 = '\0';
                if (dd - 1 >= 0 && dd - 1 <= xx) op1 = sq[dd - 1];
                if (dd - 2 >= 0 && dd - 2 <= xx) op2 = sq[dd - 2];
                if (DUPLICATE_ORIGINAL_CODE)
                {
                    if (dd - 1 < 0 && (-(dd - 1)) <= xx) op1 = sq[xx + (dd - 1) + 1];
                    if (dd - 2 < 0 && (-(dd - 2)) <= xx) op2 = sq[xx + (dd - 2) + 1];
                }
                if (dd + 1 >= 0 && dd + 1 <= xx) op3 = sq[dd + 1];
                if (dd + 2 >= 0 && dd + 2 <= xx) op4 = sq[dd + 2];
                csum += (UDF31 * (AAPARAMS[op1].UndKRH + AAPARAMS[op3].UndKRH)) + (UDF32 * (AAPARAMS[op2].UndKRH + AAPARAMS[op4].UndKRH));
            }
        }
        return csum;
    }

    private double Clusterness(string sq)
    {
        if (NOCLUSTER == 1) return 0.0;
        int bufferSize = sq.Length + 2;
        Span<char> cc = bufferSize <= 128 ? stackalloc char[bufferSize] : new char[bufferSize];
        cc[0] = '0'; cc[^1] = '0';
        for (int i = 0; i < sq.Length; i++) cc[i + 1] = EncodeClusterChar(sq[i]);
        double score = 0.0;
        foreach (var pair in CLUSTCOMB)
        {
            int occurs = CountPatternOccurrences(cc, pair.Key);
            if (occurs > 0) score += pair.Value * occurs;
        }
        return score * KSCALE;
    }

    private static char EncodeClusterChar(char aa) => aa switch
    {
        'L' or 'I' or 'W' => '5',
        'A' or 'M' or 'Y' or 'V' => '1',
        _ => '0'
    };

    private static int CountPatternOccurrences(ReadOnlySpan<char> text, string pattern)
    {
        int count = 0, index = 0;
        while (index <= text.Length - pattern.Length)
        {
            int found = text.Slice(index).IndexOf(pattern.AsSpan(), StringComparison.Ordinal);
            if (found < 0) break;
            count++; index += found + pattern.Length;
        }
        return count;
    }

    private static double Proline(string sq) =>
        sq.Contains("PPPP") ? PPPPSCORE : sq.Contains("PPP") ? PPPSCORE : sq.Contains("PP") ? PPSCORE : 0.0;

    private static double Length_scale(int sqlen) =>
        sqlen < SPLim ? 1.0 + SPSFac * (SPLim - sqlen) : sqlen > LPLim ? 1.0 / (1.0 + LPSFac * (sqlen - LPLim)) : 1.0;

    private static double Partial_charge(double pK, double pH) { double cr = Math.Pow(10.0, pK - pH); return cr / (cr + 1.0); }

    private double Electric(string sq)
    {
        int[] aaCNT = { 0, 0, 0, 0, 0, 0, 0 };
        int ss = sq.Length;
        double pk0 = AAPARAMS[sq[0]].CT, pk1 = AAPARAMS[sq[ss - 1]].NT;
        for (int i = 0; i < ss; i++) { int idx = EMap[sq[i]]; if (idx >= 0) aaCNT[idx]++; }
        double best = 0.0, min = 100000; const double step1 = 0.3;
        for (double z = 0.01; z <= 14.0; z += step1) { double check = Math.Abs(CalcR(z, pk0, pk1, aaCNT)); if (check < min) { min = check; best = z; } }
        double best1 = best; min = 100000;
        for (double z = best1 - step1; z <= best1 + step1; z += 0.01) { double check = Math.Abs(CalcR(z, pk0, pk1, aaCNT)); if (check < min) { min = check; best = z; } }
        return best;
    }

    private double CalcR(double pH, double PK0, double PK1, int[] CNTref) =>
        Partial_charge(PK0, pH)
        + CNTref[EMap['K']] * Partial_charge(AAPARAMS['K'].PK, pH)
        + CNTref[EMap['R']] * Partial_charge(AAPARAMS['R'].PK, pH)
        + CNTref[EMap['H']] * Partial_charge(AAPARAMS['H'].PK, pH)
        - CNTref[EMap['D']] * Partial_charge(pH, AAPARAMS['D'].PK)
        - CNTref[EMap['E']] * Partial_charge(pH, AAPARAMS['E'].PK)
        - CNTref[EMap['Y']] * Partial_charge(pH, AAPARAMS['Y'].PK)
        - Partial_charge(pH, PK1);

    private double NewIso(string sq, double tsum)
    {
        if (NOELECTRIC == 1)
            return 0.0;

        double mass = 0.0;
        foreach (char c in sq)
        {
            mass += AAPARAMS[c].AMASS;
        }

        double pi1 = Electric(sq);
        double lmass = 1.8014 * Math.Log(mass);
        double delta1 = pi1 - 19.107 + lmass;

        if (delta1 < 0.0)
        {
            return (tsum * Z01 + Z02) * NDELTAWT * delta1;
        }
        else if (delta1 > 0.0)
        {
            // Note: Currently returns 0 because Z03 = Z04 = 0.0
            // This is intentional per the original SSRCalc algorithm
            return (tsum * Z03 + Z04) * PDELTAWT * delta1;
        }

        return 0.0;
    }

    private static double Heli1TermAdj(string ss1, int ix2, int sqlen)
    {
        int where = 0;
        for (int i = 0; i < ss1.Length; i++) { char m = ss1[i]; if (m == 'O' || m == 'U') { where = i; if (!DUPLICATE_ORIGINAL_CODE) break; } }
        where += ix2;
        if (where < 2) return 0.20; if (where < 3) return 0.25; if (where < 4) return 0.45;
        if (where > sqlen - 3) return 0.2; if (where > sqlen - 4) return 0.75; if (where > sqlen - 5) return 0.65;
        return 1.0;
    }

    private double Helicity1(string sq)
    {
        if (NOHELIX1 == 1) return 0.0;
        int sqlen = sq.Length;
        Span<char> hcSpan = sqlen <= 128 ? stackalloc char[sqlen] : new char[sqlen];
        for (int j = 0; j < sqlen; j++) hcSpan[j] = EncodeHelicity1Char(sq[j]);
        string hc = new string(hcSpan);
        double sum = 0.0;
        for (int i = 0; i < sqlen - 3; i++)
        {
            if (i + 6 <= sqlen)
            {
                string sub6 = hc.Substring(i, 6);
                if (HlxScore6.TryGetValue(sub6, out double sc6) && sc6 > 0)
                {
                    sum += sc6 * Heli1TermAdj(sub6, i, sqlen);
                    i++;
                    continue;
                }
            }
            if (i + 5 <= sqlen)
            {
                string sub5 = hc.Substring(i, 5);
                if (HlxScore5.TryGetValue(sub5, out double sc5) && sc5 > 0)
                {
                    sum += sc5 * Heli1TermAdj(sub5, i, sqlen);
                    i++;
                    continue;
                }
            }
            if (i + 4 <= sqlen)
            {
                string sub4 = hc.Substring(i, 4);
                if (HlxScore4.TryGetValue(sub4, out double sc4) && sc4 > 0)
                {
                    sum += sc4 * Heli1TermAdj(sub4, i, sqlen);
                    i++;
                }
            }
        }
        return HELIX1SCALE * sum;
    }

    private static char EncodeHelicity1Char(char aa) => aa switch
    {
        'P' or 'H' or 'R' or 'K' => 'z',
        'W' or 'F' or 'I' or 'L' => 'X',
        'Y' or 'M' or 'V' or 'A' => 'Z',
        'D' or 'E' => 'O',
        'G' or 'S' or 'C' or 'N' or 'Q' or 'T' => 'U',
        _ => aa
    };

    private double EvalH2pattern(string pattern, string testsq, int posn, char etype)
    {
        char f01 = pattern[0];
        double prod1 = AAPARAMS[f01].H2BASCORE;
        const int OFF1 = 2; int acount = 1; char far1 = '\0', far2 = '\0';
        char testAAl = testsq[OFF1 + posn], testAAr = testsq[OFF1 + posn + 2];
        string testsqCopy = testsq.Substring(OFF1 + posn + 1);
        double mult = Connector(f01, testAAl, testAAr, "--", far1, far2);
        prod1 *= mult; if (etype == '*') prod1 *= 25.0; if (mult <= 0.0) return 0.0;
        for (int i = 1; i < pattern.Length - 2; i += 3)
        {
            string fpart = pattern.Substring(i, 2);
            char gpart = (i + 2) < pattern.Length ? pattern[i + 2] : '\0';
            double s3 = AAPARAMS[gpart].H2BASCORE;
            int iss = 0;
            if (fpart == "--") { iss = 0; far1 = '\0'; far2 = '\0'; }
            if (fpart == "<-") { iss = 1; far1 = testsqCopy[i + 1]; far2 = '\0'; }
            if (fpart == "->") { iss = -1; far1 = '\0'; far2 = testsqCopy[i + 3]; }
            testAAl = testsqCopy[i + 1 + iss]; testAAr = testsqCopy[i + 3 + iss];
            mult = Connector(gpart, testAAl, testAAr, fpart, far1, far2);
            if (etype == '*' && (mult > 0.0 || acount < 3)) prod1 = prod1 * 25.0 * s3 * mult;
            if (etype == '+') prod1 = prod1 + s3 * mult;
            if (mult <= 0.0) return prod1;
            acount++;
        }
        return prod1;
    }

    private double Connector(char acid, char lp, char rp, string ct, char far1, char far2)
    {
        double mult = 1.0;
        if (ct.Contains("<-")) mult *= 0.2; if (ct.Contains("->")) mult *= 0.1;
        mult *= AAPARAMS[lp].H2CMULT; if (lp != rp) mult *= AAPARAMS[rp].H2CMULT;
        if (acid == 'A' || acid == 'Y' || acid == 'V' || acid == 'M')
        { if (lp == 'P' || lp == 'G' || rp == 'P' || rp == 'G') mult = 0.0; if (ct.Contains("->") || ct.Contains("<-")) mult = 0.0; }
        if (acid == 'L' || acid == 'W' || acid == 'F' || acid == 'I')
        { if (((lp == 'P' || lp == 'G') || (rp == 'P' || rp == 'G')) && !ct.Contains("--")) mult = 0.0; if (((far1 == 'P' || far1 == 'G') || (far2 == 'P' || far2 == 'G')) && (ct.Contains("<-") || ct.Contains("->"))) mult = 0.0; }
        return mult;
    }

    private double[] Heli2Calc(string sq)
    {
        double[] ret = new double[2];
        if (sq.Length < 11) { ret[HISC] = 0.0; ret[GSC] = 0.0; return ret; }
        string prechop = sq, sqCopy = sq.Substring(2, sq.Length - 4);
        string pass1 = sqCopy.ReplaceAAs("WFILYMVA", "1").ReplaceAAs("GSPCNKQHRTDE", "0");
        string best = ""; double hiscore = 0.0; int best_pos = 0;
        for (int i = 0; i < pass1.Length; i++)
        {
            if (pass1[i] == '1')
            {
                string lc = pass1.Substring(i), sq2 = sqCopy.Substring(i);
                StringBuilder patBuilder = new();
                int zap = 0, subt = 0;
                while (zap <= 50 && subt < 2)
                {
                    char f1 = zap >= 0 && zap < lc.Length ? lc[zap] : '0';
                    char f2 = zap - 1 >= 0 && zap - 1 < lc.Length ? lc[zap - 1] : '0';
                    char f3 = zap + 1 >= 0 && zap + 1 < lc.Length ? lc[zap + 1] : '0';
                    if (f1 == '1') { if (zap > 0) patBuilder.Append("--"); patBuilder.Append(sq2[zap]); }
                    else if (f2 == '1' && f1 == '0') { subt++; if (subt < 2) { patBuilder.Append("->"); patBuilder.Append(sq2[zap - 1]); } }
                    else if (f3 == '1' && f1 == '0') { subt++; if (subt < 2) { patBuilder.Append("<-"); patBuilder.Append(sq2[zap + 1]); } }
                    if (f1 == '0' && f2 == '0' && f3 == '0') zap = 1000;
                    zap += 3;
                }
                if (patBuilder.Length > 4)
                {
                    string pat = patBuilder.ToString();
                    double skore = EvalH2pattern(pat, prechop, i - 1, '*');
                    if (skore >= hiscore) { hiscore = skore; best = pat; best_pos = i; }
                }
            }
        }
        if (hiscore > 0.0)
        {
            double gscore = hiscore;  // Preserve '*' score for GSC (used in forward/backward comparison)
            hiscore = EvalH2pattern(best, prechop, best_pos - 1, '+');  // '+' score for HISC
            ret[HISC] = hiscore;
            ret[GSC] = gscore;
            return ret;
        }
        ret[HISC] = 0.0; ret[GSC] = 0.0; return ret;
    }

    private double Helicity2(string sq)
    {
        if (NOHELIX2 == 1) return 0.0;
        string Bksq = sq.Backwards();
        double[] fhg = Heli2Calc(sq), rhg = Heli2Calc(Bksq);
        double h2FwBk = rhg[GSC] > fhg[GSC] ? rhg[HISC] : fhg[HISC];
        double lenMult = sq.Length > 30 ? 1 : 0, NoPMult = sq.Contains('P') ? 0.0 : 0.75;
        return HELIX2SCALE * (1.0 + lenMult + NoPMult) * h2FwBk;
    }

    private double Helectric(string sq)
    {
        if (NOEHEL == 1 || sq.Length > 14 || sq.Length < 4) return 0.0;
        string mpart = sq.Substring(sq.Length - 4);
        if (mpart[0] == 'D' || mpart[0] == 'E')
        {
            mpart = mpart.Substring(1, 2);
            if (mpart.ContainsAA("PGKRH")) return 0.0;
            mpart = mpart.ReplaceAAs("LI", "X").ReplaceAAs("AVYFWM", "Z").ReplaceAAs("GSPCNKQHRTDE", "U");
            return mpart switch { "XX" => 1.0, "ZX" => 0.5, "XZ" => 0.5, "ZZ" => 0.4, "XU" => 0.4, "UX" => 0.4, "ZU" => 0.2, "UZ" => 0.2, _ => 0 };
        }
        return 0;
    }

    public class AAParams
    {
        public double RC { get; private set; }
        public double RC1 { get; private set; }
        public double RC2 { get; private set; }
        public double RCN { get; private set; }
        public double RCN2 { get; private set; }
        public double RCS { get; private set; }
        public double RC1S { get; private set; }
        public double RC2S { get; private set; }
        public double RCNS { get; private set; }
        public double RCN2S { get; private set; }
        public double UndKRH { get; private set; }
        public double AMASS { get; private set; }
        public double CT { get; private set; }
        public double NT { get; private set; }
        public double PK { get; private set; }
        public double H2BASCORE { get; private set; }
        public double H2CMULT { get; private set; }
        public AAParams(double rc, double rc1, double rc2, double rcn, double rcn2, double rcs, double rc1s, double rc2s, double rcns, double rcn2s, double undkrh, double amass, double ct, double nt, double pk, double h2bascore, double h2cmult)
        { RC = rc; RC1 = rc1; RC2 = rc2; RCN = rcn; RCN2 = rcn2; RCS = rcs; RC1S = rc1s; RC2S = rc2s; RCNS = rcns; RCN2S = rcn2s; UndKRH = undkrh; AMASS = amass; CT = ct; NT = nt; PK = pk; H2BASCORE = h2bascore; H2CMULT = h2cmult; }
    }
}

internal static class HelpersLocal
{
    public static string ReplaceAAs(this IEnumerable<char> s, string aas, string newValue)
    {
        StringBuilder sb = new();
        bool allAAs = aas == "A-Z";

        foreach (char c in s)
        {
            if (!allAAs && aas.IndexOf(c) != -1)
            {
                sb.Append(newValue);
            }
            else if (allAAs && char.IsLetter(c) && char.IsUpper(c))
            {
                sb.Append(newValue);
            }
            else
            {
                sb.Append(c);
            }
        }

        return sb.ToString();
    }
    /// <summary>
    /// Inspects a sequence of amino acids, and returns true if it contains
    /// any of the designated amino acid characters.
    /// </summary>
    /// <param name="s">Amino acid sequence</param>
    /// <param name="aas">List of characters to search for</param>
    /// <returns>True if any of the amino acid characters are found</returns>
    public static bool ContainsAA(this IEnumerable<char> s, string aas)
    {
        ReadOnlySpan<char> aasSpan = aas.AsSpan();
        foreach (char c in s)
        {
            if (aasSpan.Contains(c))
                return true;
        }
        return false;
    }

    public static string Backwards(this IEnumerable<char> s)
    {
        StringBuilder sb = new();
        foreach (char c in s.Reverse())
        {
            sb.Append(c);
        }
        return sb.ToString();
    }
}