using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TestUniDec
{
	public partial class UniDecClasses
	{
		public struct Config
		{
			string infile;
			string outfile;
			int numit;
			int numz;
			int endz;
			int startz;
			float zsig;
			float psig;
			float beta;
			float mzsig;
			float msig;
			float molig;
			float massub;
			float masslb;
			int psfun;
			float mtabsig;
			string mfile;
			string manualfile;
			int mflag;
			float massbins;
			int limitflag;
			float psthresh;
			int speedyflag;
			int linflag;
			int aggressiveflag;
			float adductmass;
			int rawflag;
			float nativezub;
			float nativezlb;
			int poolflag;
			int manualflag;
			float intthresh;
			float peakshapeinflate;
			float killmass;
			int fixedmassaxis;
			int isotopemode;
			int filetype;
			int imflag;
			//IM Parameters
			float dtsig;
			float csig;
			float ccsub;
			float ccslb;
			float ccsbins;
			float temp;
			float press;
			float volt;
			float tcal1;
			float tcal2;
			float tcal3;
			float tcal4;
			float twaveflag;
			float hmass;
			float to;
			float len;
			float edc;
			float nativeccsub;
			float nativeccslb;
			int baselineflag;
			int noiseflag;
			int zout;
			int metamode;
			float minmz;
			float maxmz;
			int mzbins;
			float bsub;
			float datareduction;
			float peakwin;
			float peakthresh;
			float exwindow;
			int exchoice;
			int exchoicez;
			float exthresh;
			int exnorm;
			int exnormz;
			int peaknorm;
			int orbimode;
			int datanorm;
			//Experimental Parameters
			int filterwidth;
			float zerolog;
			int lengthmz;
			int mfilelen;
			int isolength;
			// DoubleDec Parameters
			int doubledec;
			char kernel;
			int file_id;
			string dataset;
		}
	}
}
