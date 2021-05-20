using System;
using System.Collections.Generic;
using System.IO;
using System.Reflection;
using System.Text;

namespace MzLibUtil
{
    public static class TextFileReading
    {
        readonly static FieldInfo charPosField = typeof(StreamReader).GetField("_charPos", BindingFlags.NonPublic | BindingFlags.Instance | BindingFlags.DeclaredOnly);
        readonly static FieldInfo charLenField = typeof(StreamReader).GetField("_charLen", BindingFlags.NonPublic | BindingFlags.Instance | BindingFlags.DeclaredOnly);
        readonly static FieldInfo charBufferField = typeof(StreamReader).GetField("_charBuffer", BindingFlags.NonPublic | BindingFlags.Instance | BindingFlags.DeclaredOnly);

        /// <summary>
        /// The StreamReader object does not have a reliable method to get a byte position in
        /// a file. This method calculates it.
        /// 
        /// Code from: https://stackoverflow.com/questions/10189270/tracking-the-position-of-the-line-of-a-streamreader/34744440
        /// 
        /// See also alternative method:
        /// Code from: https://stackoverflow.com/questions/5404267/streamreader-and-seeking
        /// </summary>
        public static long GetByteOffsetAtCurrentPosition(StreamReader reader)
        {
            if (reader == null)
            {
                throw new MzLibException("Attempted to get position in the text stream but the StreamReader was null");
            }

            var charBuffer = (char[])charBufferField.GetValue(reader);
            var charLen = (int)charLenField.GetValue(reader);
            var charPos = (int)charPosField.GetValue(reader);

            return reader.BaseStream.Position - reader.CurrentEncoding.GetByteCount(charBuffer, charPos, charLen - charPos);
        }
    }
}
