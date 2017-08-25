// Copyright 2016 Stefan Solntsev
//
// This file (PeriodicTable.cs) is part of UsefulProteomicsDatabases.
//
// UsefulProteomicsDatabases is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// UsefulProteomicsDatabases is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with UsefulProteomicsDatabases. If not, see <http://www.gnu.org/licenses/>.

using Chemistry;
using System;
using System.Globalization;
using System.IO;
using System.Text.RegularExpressions;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// The Periodic Table of Elements.
    /// </summary>
    public static class PeriodicTableLoader
    {
        #region Public Methods

        public static void Load(string elementLocation)
        {
            using (StreamReader sr = new StreamReader(elementLocation))
            {
                string line = sr.ReadLine();
                while (!line.Contains("Atomic Number"))
                {
                    line = sr.ReadLine();
                }
                var prevAtomicNumber = -1;
                Element element = new Element("fake", prevAtomicNumber, -1);
                do
                {
                    int atomicNumber = Convert.ToInt32(Regex.Match(line, @"\d+").Value, CultureInfo.InvariantCulture);

                    line = sr.ReadLine();
                    string atomicSymbol = Regex.Match(line, @"[A-Za-z]+$").Value;

                    line = sr.ReadLine();
                    int massNumber = Convert.ToInt32(Regex.Match(line, @"\d+").Value, CultureInfo.InvariantCulture);

                    line = sr.ReadLine();
                    double atomicMass = Convert.ToDouble(Regex.Match(line, @"[\d\.]+").Value, CultureInfo.InvariantCulture);

                    line = sr.ReadLine();
                    double abundance;
                    if (Regex.Match(line, @"[\d\.]+").Success)
                    {
                        abundance = Convert.ToDouble(Regex.Match(line, @"[\d\.]+").Value, CultureInfo.InvariantCulture);
                    }
                    else
                    {
                        sr.ReadLine();
                        sr.ReadLine();
                        sr.ReadLine();
                        line = sr.ReadLine();
                        continue;
                    }

                    line = sr.ReadLine();
                    double averageMass;
                    if (Regex.Match(line, @"\[").Success)
                    {
                        double averageMass1 = Convert.ToDouble(Regex.Match(line, @"(?<=\[)[\d\.]+").Value, CultureInfo.InvariantCulture);
                        var kkajsdf = Regex.Match(line, @"(?<=,)[\d\.]+").Value;
                        double averageMass2 = Convert.ToDouble(kkajsdf, CultureInfo.InvariantCulture);
                        averageMass = (averageMass1 + averageMass2) / 2;
                    }
                    else
                        averageMass = Convert.ToDouble(Regex.Match(line, @"[\d\.]+").Value, CultureInfo.InvariantCulture);

                    if (atomicNumber != prevAtomicNumber)
                    {
                        element = new Element(atomicSymbol, atomicNumber, averageMass);
                        PeriodicTable.Add(element);
                    }

                    element.AddIsotope(massNumber, atomicMass, abundance);

                    sr.ReadLine();
                    sr.ReadLine();
                    line = sr.ReadLine();
                    prevAtomicNumber = atomicNumber;
                } while (line.Contains("Atomic Number"));
            }
        }

        #endregion Public Methods
    }
}