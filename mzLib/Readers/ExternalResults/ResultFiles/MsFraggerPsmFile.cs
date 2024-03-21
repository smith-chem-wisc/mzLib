using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers
{
    public class MsFraggerPsmFile : ResultFile<MsFraggerPsm>, IResultFile
    {
        public override SupportedFileType FileType => SupportedFileType.MsFraggerPsm;
        public override Software Software { get; set; }
        public MsFraggerPsmFile(string filePath) : base(filePath, Software.MsFragger) { }

        public override void LoadResults()
        {
            throw new NotImplementedException();
        }

        public override void WriteResults(string outputPath)
        {
            throw new NotImplementedException();
        }
    }
}
