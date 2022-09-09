using System.Linq;

namespace MzLibUtil
{
    public class PeriodTolerantFilenameWithoutExtension
    {
        public static string GetPeriodTolerantFilenameWithoutExtension(string filepathAndOrName)
        {
            string windowsSplit = filepathAndOrName.Split('\\').Last();
            string linuxSplit = filepathAndOrName.Split('/').Last();

            if(windowsSplit.Length <= linuxSplit.Length)
            {
                int periodBeforeExtension = windowsSplit.LastIndexOf('.');
                if (periodBeforeExtension > 0)
                {
                    return windowsSplit[..periodBeforeExtension];
                }
                else
                {
                    return windowsSplit;
                }  
            }
            else
            {
                int periodBeforeExtension = linuxSplit.LastIndexOf('.');
                if (periodBeforeExtension > 0)
                {
                    return linuxSplit[..periodBeforeExtension];
                }
                else
                {
                    return linuxSplit;
                }
            }
        }
    }
}
