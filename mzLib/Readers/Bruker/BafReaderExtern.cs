using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System;
using System.Runtime.InteropServices;


namespace Readers
{
    /// <summary>
    /// Methods in this class are all from Bruker. 
    /// </summary>
    internal class BafReaderExtern
    {
        [DllImport("baf2sql_c", CallingConvention = CallingConvention.Cdecl, EntryPoint = "baf2sql_get_sqlite_cache_filename")]
        public static extern UInt32 baf2sql_get_sqlite_cache_filename
            (byte[] sql_filename_buf_utf8, UInt32 sql_filename_buflen,
                byte[] baf_filename_utf8);

        [DllImport("baf2sql_c", CallingConvention = CallingConvention.Cdecl)]
        public static extern UInt64 baf2sql_array_open_storage
               (int ignore_calibrator_ami, byte[] filename_utf8);

        [DllImport("baf2sql_c", CallingConvention = CallingConvention.Cdecl)]
        public static extern void baf2sql_array_close_storage
               (UInt64 handle);

        [DllImport("baf2sql_c", CallingConvention = CallingConvention.Cdecl)]
        public static extern void baf2sql_array_get_num_elements
               (UInt64 handle, UInt64 id, ref UInt64 num_elements);

        [DllImport("baf2sql_c", CallingConvention = CallingConvention.Cdecl)]
        public static extern int baf2sql_array_read_double
               (UInt64 handle, UInt64 id, double[] buf);

        [DllImport("baf2sql_c", CallingConvention = CallingConvention.Cdecl)]
        public static extern int baf2sql_array_read_float
               (UInt64 handle, UInt64 id, float[] buf);

        [DllImport("baf2sql_c", CallingConvention = CallingConvention.Cdecl)]
        public static extern int baf2sql_array_read_uint32
               (UInt64 handle, UInt64 id, UInt32[] buf);

        [DllImport("baf2sql_c", CallingConvention = CallingConvention.Cdecl)]
        public static extern UInt32 baf2sql_get_last_error_string
               (StringBuilder buf, UInt32 len);

        [DllImport("baf2sql_c", CallingConvention = CallingConvention.Cdecl)]
        public static extern void baf2sql_set_num_threads
               (UInt32 n);

        /* ----------------------------------------------------------------------------------------------- */

        [Serializable()]
        public class Baf2SqlException : System.Exception
        {
            public Baf2SqlException() : base() { }
            public Baf2SqlException(string message) : base(message) { }
            public Baf2SqlException(string message, System.Exception inner) : base(message, inner) { }
            protected Baf2SqlException(System.Runtime.Serialization.SerializationInfo info,
                System.Runtime.Serialization.StreamingContext context)
            { }
        }

        /* Throw last error string as an exception. */
        public static void ThrowLastBaf2SqlError()
        {
            StringBuilder buf = new StringBuilder("");
            UInt32 len = baf2sql_get_last_error_string(buf, 0);
            buf.EnsureCapacity((int)(len + 1));
            baf2sql_get_last_error_string(buf, len);
            throw new Baf2SqlException(buf.ToString());
        }

        /* ----------------------------------------------------------------------------------------------- */

        /* Find out the file name of the SQL cache corresponding to the specified BAF file.
         * (If the SQL cache doesn't exist yet, it will be created.) */
        public static String GetSQLiteCacheFilename(String baf_filename)
        {
            byte[] buf = new byte[0];
            byte[] baf_filename_utf8 = ConvertStringToUTF8ByteArray(baf_filename);

            UInt32 len = baf2sql_get_sqlite_cache_filename(buf, 0, baf_filename_utf8);
            if (len == 0) ThrowLastBaf2SqlError();

            buf = new byte[len];
            len = baf2sql_get_sqlite_cache_filename(buf, len, baf_filename_utf8);
            if (len == 0) ThrowLastBaf2SqlError();

            return Encoding.UTF8.GetString(buf, 0, buf.Length - 1);
        }

        /* ----------------------------------------------------------------------------------------------- */

        /* Given the Id of one spectral component (e.g., a 'ProfileMzId' from the SQL cache),
         * load the binary data from the BAF (returning a double array). */
        public static double[] GetBafDoubleArray(UInt64 handle, UInt64 id)
        {
            UInt64 n = 0;
            baf2sql_array_get_num_elements(handle, id, ref n);

            double[] myArray = new double[n];
            int rc = baf2sql_array_read_double(handle, id, myArray);
            if (rc == 0) ThrowLastBaf2SqlError();

            return myArray;
        }

        /* Return array 'id', converting to float format */
        public static float[] GetBafFloatArray(UInt64 handle, UInt64 id)
        {
            UInt64 n = 0;
            baf2sql_array_get_num_elements(handle, id, ref n);

            float[] myArray = new float[n];
            int rc = baf2sql_array_read_float(handle, id, myArray);
            if (rc == 0) ThrowLastBaf2SqlError();

            return myArray;
        }

        /* Return array 'id', converting to UInt32 format */
        public static UInt32[] GetBafUInt32Array(UInt64 handle, UInt64 id)
        {
            UInt64 n = 0;
            baf2sql_array_get_num_elements(handle, id, ref n);

            UInt32[] myArray = new UInt32[n];
            int rc = baf2sql_array_read_uint32(handle, id, myArray);
            if (rc == 0) ThrowLastBaf2SqlError();

            return myArray;
        }
        public static byte[] ConvertStringToUTF8ByteArray(String input)
        {
            byte[] utf8 = Encoding.UTF8.GetBytes(input);
            var result = new byte[utf8.Length];
            Array.Copy(utf8, result, utf8.Length);

            return result;
        }
    }
}
