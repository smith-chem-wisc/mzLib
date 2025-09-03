using MzLibUtil;
using System.Data.SQLite;

namespace Readers
{
    internal enum TimsTofMsMsType
    {
        MS = 0,
        MSMSFragment = 2,
        PASEF = 8,
        DIA = 9,
        PRM = 10
    }

    internal enum TimsTofAcquisitionMode
    {
        MS = 0,
        AutoMSMS = 1,
        MRM = 2,
        inSourceCID = 3,
        broadbandCID = 4,
        PASEF = 8,
        DIA = 9,
        PRM = 10,
        Maldi = 20
    }

    /// <summary>
    /// This class stores information take from the .tdf SQLite database file
    /// Every frame in the file has 9 pieces of metadata that can be accessed by 
    /// selecting the appropriate array. All arrays are zero-based!!!
    /// EX: ScanMode[0] will return the scan mode of the first frame (FrameID = 1) in the file
    /// </summary>
    internal class FrameTable
    {
        internal long[] OneBasedFrameIndex { get; private set; }
        internal char[] Polarity { get; private set; }
        internal int[] NumScans { get; private set; }
        internal int[] ScanMode { get; private set; }
        internal int[] MsMsType { get; private set; }
        internal int[] TotalNumberOfPeaks { get; private set; }
        internal int[] TotalIntensity { get; private set; }
        internal float[] RetentionTime { get; private set; }
        internal float[] FillTime { get; private set; }

        internal TimsTofMsMsType GetAnalysisType(int frameId)
        {
            if (frameId == 0 || frameId > MsMsType.Length) throw new IndexOutOfRangeException("Invalid frame ID!");
            if (MsMsType[frameId - 1].ToEnum<TimsTofMsMsType>(out var analysisType))
                return analysisType;
            else
                throw new MzLibException("Unrecognized MS/MS method.");
        }

        internal FrameTable(SQLiteConnection connection, int numberOfRows, TimsTofFileType fileType = TimsTofFileType.TDF)
        {
            switch (fileType)
            {
                case TimsTofFileType.TDF:
                    PopulateTableForTdf(connection, numberOfRows);
                    break;
                case TimsTofFileType.TSF:
                    PopulateTableForTsf(connection, numberOfRows);
                    break;
            }
        }

        private void PopulateTableForTdf(SQLiteConnection connection, int numberOfRows)
        {
            using var command = new SQLiteCommand(connection);
            command.CommandText = @"SELECT f.Id, f.Polarity, f.NumScans," +
                                  " f.ScanMode, f.MsMsType, f.NumPeaks, f.SummedIntensities," +
                                  " f.Time, f.AccumulationTime FROM Frames f;";
            using var reader = command.ExecuteReader();

            OneBasedFrameIndex = new long[numberOfRows];
            Polarity = new char[numberOfRows];
            NumScans = new int[numberOfRows];
            ScanMode = new int[numberOfRows];
            MsMsType = new int[numberOfRows];
            TotalNumberOfPeaks = new int[numberOfRows];
            TotalIntensity = new int[numberOfRows];
            RetentionTime = new float[numberOfRows];
            FillTime = new float[numberOfRows];

            // Populate arrays by reading in the  table
            for (int i = 0; i < numberOfRows; i++)
            {
                if (!reader.Read()) break;
                OneBasedFrameIndex[i] = reader.GetInt64(0);
                Polarity[i] = reader.GetString(1)[0];
                NumScans[i] = reader.GetInt32(2);
                ScanMode[i] = reader.GetInt32(3);
                MsMsType[i] = reader.GetInt32(4);
                TotalNumberOfPeaks[i] = reader.GetInt32(5);
                TotalIntensity[i] = reader.GetInt32(6);
                RetentionTime[i] = reader.GetFloat(7);
                FillTime[i] = reader.GetFloat(8);
            }
        }

        private void PopulateTableForTsf(SQLiteConnection connection, int numberOfRows)
        {
            using var command = new SQLiteCommand(connection);
            command.CommandText = @"SELECT f.Id, f.Polarity," +
                                  " f.ScanMode, f.MsMsType, f.NumPeaks, f.SummedIntensities," +
                                  " f.Time FROM Frames f;";
            using var reader = command.ExecuteReader();

            OneBasedFrameIndex = new long[numberOfRows];
            Polarity = new char[numberOfRows];
            NumScans = null;
            ScanMode = new int[numberOfRows];
            MsMsType = new int[numberOfRows];
            TotalNumberOfPeaks = new int[numberOfRows];
            TotalIntensity = new int[numberOfRows];
            RetentionTime = new float[numberOfRows];
            FillTime = null;

            // Populate arrays by reading in the  table
            for (int i = 0; i < numberOfRows; i++)
            {
                if (!reader.Read()) break;
                OneBasedFrameIndex[i] = reader.GetInt64(0);
                Polarity[i] = reader.GetString(1)[0];
                ScanMode[i] = reader.GetInt32(2);
                MsMsType[i] = reader.GetInt32(3);
                TotalNumberOfPeaks[i] = reader.GetInt32(4);
                TotalIntensity[i] = reader.GetInt32(5);
                RetentionTime[i] = reader.GetFloat(6);
            }
        }
    }
}
