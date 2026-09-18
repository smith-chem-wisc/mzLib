using System;
using System.Buffers;
using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;
using Quantification.Interfaces;

namespace Quantification.Strategies
{
    /// <summary>
    /// Shared machinery for the collapse strategies. A subclass only declares which of the three
    /// sample dimensions it collapses; grouping, aggregation and the identity of the surviving column
    /// are handled here.
    ///
    /// Columns are grouped by <see cref="ISampleInfo.Condition"/> plus whichever of biological
    /// replicate, fraction and technical replicate the subclass does *not* collapse. Groups appear in
    /// the order their first member appeared, so output column order is deterministic.
    /// </summary>
    public abstract class CollapseStrategyBase : ICollapseStrategy
    {
        public abstract string Name { get; }

        protected abstract bool CollapsesBiologicalReplicate { get; }
        protected abstract bool CollapsesFraction { get; }
        protected abstract bool CollapsesTechnicalReplicate { get; }

        public QuantMatrix<T> CollapseSamples<T>(QuantMatrix<T> quantMatrix, IAggregationStrategy aggregation)
            where T : IEquatable<T>
        {
            if (aggregation == null) throw new ArgumentNullException(nameof(aggregation));

            var groups = GroupColumns(quantMatrix.ColumnKeys);

            var collapsedMatrix = new QuantMatrix<T>(
                quantMatrix.RowKeys,
                groups.Select(g => CollapsedSampleInfo(g[0].Column)).ToList(),
                quantMatrix.ExperimentalDesign);

            int widestGroup = groups.Count == 0 ? 1 : groups.Max(g => g.Count);
            double[] collapsedRow = ArrayPool<double>.Shared.Rent(Math.Max(groups.Count, 1));
            double[] groupValues = ArrayPool<double>.Shared.Rent(widestGroup);

            try
            {
                for (int row = 0; row < quantMatrix.RowCount; row++)
                {
                    for (int groupIndex = 0; groupIndex < groups.Count; groupIndex++)
                    {
                        var group = groups[groupIndex];
                        for (int i = 0; i < group.Count; i++)
                            groupValues[i] = quantMatrix.Matrix[row, group[i].Index];

                        collapsedRow[groupIndex] = aggregation.Aggregate(groupValues.AsSpan(0, group.Count));
                    }

                    collapsedMatrix.SetRow(quantMatrix.RowKeys[row], collapsedRow);
                }
            }
            finally
            {
                ArrayPool<double>.Shared.Return(collapsedRow, clearArray: true);
                ArrayPool<double>.Shared.Return(groupValues, clearArray: true);
            }

            return collapsedMatrix;
        }

        /// <summary>
        /// Buckets the columns by the dimensions this strategy keeps, preserving first-appearance order
        /// both between groups and within each group.
        /// </summary>
        private List<List<(ISampleInfo Column, int Index)>> GroupColumns(IReadOnlyList<ISampleInfo> columns)
        {
            var order = new List<(string, int, int, int)>();
            var groups = new Dictionary<(string, int, int, int), List<(ISampleInfo, int)>>();

            for (int index = 0; index < columns.Count; index++)
            {
                var column = columns[index];
                var key = (
                    column.Condition ?? string.Empty,
                    CollapsesBiologicalReplicate ? 0 : column.BiologicalReplicate,
                    CollapsesFraction ? 0 : column.Fraction,
                    CollapsesTechnicalReplicate ? 0 : column.TechnicalReplicate);

                if (!groups.TryGetValue(key, out var members))
                {
                    members = new List<(ISampleInfo, int)>();
                    groups[key] = members;
                    order.Add(key);
                }
                members.Add((column, index));
            }

            return order.Select(k => groups[k]).ToList();
        }

        /// <summary>
        /// The <see cref="ISampleInfo"/> standing for a collapsed group: the dimensions this strategy
        /// collapsed are zeroed, and everything else is taken from the group's first member.
        /// </summary>
        /// <remarks>
        /// An <see cref="IsobaricQuantSampleInfo"/> collapses to another
        /// <see cref="IsobaricQuantSampleInfo"/> so that the channel label, plex, reporter m/z and
        /// reference-channel flag survive. Returning a plain <see cref="SpectraFileInfo"/> would
        /// silently disable <see cref="ReferenceChannelNormalization"/> for everything downstream.
        /// </remarks>
        internal ISampleInfo CollapsedSampleInfo(ISampleInfo source)
        {
            int biorep = CollapsesBiologicalReplicate ? 0 : source.BiologicalReplicate;
            int fraction = CollapsesFraction ? 0 : source.Fraction;
            int techrep = CollapsesTechnicalReplicate ? 0 : source.TechnicalReplicate;

            string label = $"{source.Condition}_Bio{biorep}_Frac{fraction}_Tech{techrep}_Collapsed";

            if (source is IsobaricQuantSampleInfo isobaric)
            {
                return new IsobaricQuantSampleInfo(
                    fullFilePathWithExtension: label,
                    condition: isobaric.Condition,
                    biologicalReplicate: biorep,
                    technicalReplicate: techrep,
                    fraction: fraction,
                    plexId: isobaric.PlexId,
                    channelLabel: isobaric.ChannelLabel,
                    reporterIonMz: isobaric.ReporterIonMz,
                    isReferenceChannel: isobaric.IsReferenceChannel);
            }

            return new SpectraFileInfo(
                fullFilePathWithExtension: label,
                condition: source.Condition,
                biorep: biorep,
                techrep: techrep,
                fraction: fraction);
        }
    }
}
