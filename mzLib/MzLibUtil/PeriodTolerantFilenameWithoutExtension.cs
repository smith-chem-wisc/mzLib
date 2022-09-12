using System.IO;
using System.Linq;

namespace MzLibUtil
{
    public class PeriodTolerantFilenameWithoutExtension
    {
        public static string GetPeriodTolerantFilenameWithoutExtension(string filepathAndOrName)
        {
            string filenameWithOrWithoutExtension = Path.GetFileName(filepathAndOrName);

            int periodBeforeExtension = filenameWithOrWithoutExtension.LastIndexOf('.');
            if (periodBeforeExtension > 0)
            {
                string compressedQuestionMark = filenameWithOrWithoutExtension[periodBeforeExtension..];
                if(compressedQuestionMark == ".gz" || compressedQuestionMark == ".zip")
                {
                    periodBeforeExtension = filenameWithOrWithoutExtension[..periodBeforeExtension].LastIndexOf('.');
                }  
                return filenameWithOrWithoutExtension[..periodBeforeExtension];
            }
            else
            {
                return filenameWithOrWithoutExtension;
            }  
        }
    }
}
