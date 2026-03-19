using Chemistry;
using MzLibUtil;
using Omics.Fragmentation;
using System.Globalization;
using System.Text;
using System.Text.RegularExpressions;
using Easy.Common.Extensions;
using Omics.Fragmentation.Peptide;
using Omics.SpectrumMatch;

namespace Readers.SpectralLibrary
{
	public class SpectralLibrary : ResultFile<LibrarySpectrum>, IResultFile
	{
		public override SupportedFileType FileType => FilePath.ParseFileType();
		public override Software Software { get; set; }
		public SpectralLibrary() : base() { }
		public SpectralLibrary(string filePath) : base(filePath, Software.MetaMorpheus) { }

		public override void LoadResults()
		{
			Results = GetAllLibrarySpectra().ToList();
		}

		//This is from WriteSpectrumLibrary in MetaMorpheusTask
		public override void WriteResults(string outputPath)
		{
			using (StreamWriter output = new StreamWriter(outputPath))
			{
				foreach (var x in Results)
				{
					output.WriteLine(x.ToString());
				}
			}
		}

		/// <summary>Ordered list of all library file paths supplied by the caller at construction time.</summary>
		private List<string> LibraryPaths;

		/// <summary>
		/// Maps the "Sequence/Charge" lookup key to the file path and byte offset of the corresponding
		/// spectrum record inside an MSP/pDeep/ms2pip text library file.
		/// Not used for .msl files — those are looked up directly via <see cref="_mslLibraries"/>.
		/// </summary>
		private Dictionary<string, (string filePath, long byteOffset)> SequenceToFileAndLocation;

		/// <summary>
		/// FIFO queue of lookup keys currently held in <see cref="LibrarySpectrumBuffer"/>, used to
		/// evict the oldest entry when the buffer exceeds <see cref="MaxElementsInBuffer"/>.
		/// </summary>
		private Queue<string> LibrarySpectrumBufferList;

		/// <summary>
		/// LRU cache of recently accessed <see cref="LibrarySpectrum"/> objects, keyed by
		/// "Sequence/Charge".  Avoids repeated disk seeks for frequently queried entries.
		/// </summary>
		public Dictionary<string, LibrarySpectrum> LibrarySpectrumBuffer;

		/// <summary>Maximum number of spectra held in <see cref="LibrarySpectrumBuffer"/> at once.</summary>
		private int MaxElementsInBuffer = 10000;

		/// <summary>
		/// One open <see cref="StreamReader"/> per text library file path, used for random-access
		/// reads via byte-offset seeks.  Not used for .msl files.
		/// </summary>
		private Dictionary<string, StreamReader> StreamReaders;

		/// <summary>
		/// One loaded <see cref="MslLibrary"/> per .msl file path.
		/// Populated in the constructor for every path whose extension is <c>.msl</c>.
		/// Disposed in <see cref="CloseConnections"/>.
		/// </summary>
		private Dictionary<string, MslLibrary> _mslLibraries;

		/// <summary>
		/// Regex for standard terminal fragment ion annotations, e.g. "b5", "y12^2", "y8-17.0265".
		/// Group 1 = ion type letters (e.g. "b", "y", "zDot").
		/// Group 2 = fragment series number (digits).
		/// Group 3 = charge (digits after "^") or neutral-loss sign + digits, or empty.
		/// Does NOT match internal ion annotations — those are handled by <see cref="InternalIonRegex"/>.
		/// </summary>
		private static Regex IonParserRegex = new Regex(@"^(\D{1,})(\d{1,})(?:[\^]|$)(-?\d{1,}|$)");

		/// <summary>
		/// Regex for internal fragment ion annotations produced by
		/// <see cref="Omics.SpectralMatch.MslSpectralLibrary.MslFragmentIon.Annotation"/>.
		/// Format: <c>{PrimaryType}I{SecondaryType}[{start}-{end}]</c> with an optional
		/// <c>^{charge}</c> suffix, e.g. <c>bIb[3-6]</c> or <c>aIb[2-5]^2</c>.
		/// Groups:
		///   1 = primary product type string (e.g. "b", "a")
		///   2 = secondary product type string (e.g. "b", "y")
		///   3 = start residue number
		///   4 = end residue number
		///   5 = charge (optional; absent means charge 1)
		/// </summary>
		private static Regex InternalIonRegex =
			new Regex(@"^([a-zA-Z]+)I([a-zA-Z]+)\[(\d+)-(\d+)\](?:\^(\d+))?");

		private static Dictionary<string, string> PrositToMetaMorpheusModDictionary = new Dictionary<string, string>
		{
			{ "Oxidation","[Common Variable:Oxidation on M]" },
			{ "Carbamidomethyl", "[Common Fixed:Carbamidomethyl on C]" }
		};

		private static Dictionary<string, string> pDeepToMetaMorpheusModDictionary = new Dictionary<string, string>
		{
			{ "Oxidation","[Common Variable:Oxidation on M]" },
			{"CAM", "[Common Fixed:Carbamidomethyl on C]" }
		};

		/// <summary>
		/// Opens all libraries in <paramref name="pathsToLibraries"/>.
		/// .msl files are opened via <see cref="MslFileTypeHandler.Open"/> and stored in
		/// <see cref="_mslLibraries"/>; all other formats are indexed into
		/// <see cref="SequenceToFileAndLocation"/> via <see cref="IndexSpectralLibrary"/>.
		/// </summary>
		/// <param name="pathsToLibraries">
		///   One or more absolute or relative paths to spectral library files.
		///   May be a mix of .msl and text-based (.msp, pdeep, ms2pip) files.
		/// </param>
		public SpectralLibrary(List<string> pathsToLibraries)
		{
			LibraryPaths = pathsToLibraries;
			SequenceToFileAndLocation = new Dictionary<string, (string, long)>();
			LibrarySpectrumBufferList = new Queue<string>();
			LibrarySpectrumBuffer = new Dictionary<string, LibrarySpectrum>();
			StreamReaders = new Dictionary<string, StreamReader>();
			_mslLibraries = new Dictionary<string, MslLibrary>();

			foreach (var path in LibraryPaths)
			{
				if (MslFileTypeHandler.IsMslFile(path))
				{
					// .msl files are binary; open via the MSL reader rather than a StreamReader
					_mslLibraries[path] = MslFileTypeHandler.Open(path);
				}
				else
				{
					// Text-based formats: build the byte-offset index as before
					IndexSpectralLibrary(path);
				}
			}
		}

		/// <summary>
		/// Returns true when any loaded library contains a spectrum for the given
		/// <paramref name="sequence"/> and <paramref name="charge"/>.
		/// Checks .msl libraries first, then the text-library index.
		/// </summary>
		/// <param name="sequence">Modified sequence in mzLib bracket notation.</param>
		/// <param name="charge">Precursor charge state.</param>
		/// <returns>True if the spectrum is present in any loaded library; false otherwise.</returns>
		public bool ContainsSpectrum(string sequence, int charge)
		{
			// Check MSL libraries first
			foreach (var mslLib in _mslLibraries.Values)
			{
				if (mslLib.TryGetLibrarySpectrum(sequence, charge, out _))
					return true;
			}

			string lookupString = sequence + "/" + charge;
			return SequenceToFileAndLocation.ContainsKey(lookupString);
		}

		/// <summary>
		/// Retrieves the <see cref="LibrarySpectrum"/> for the given <paramref name="sequence"/>
		/// and <paramref name="charge"/>, searching .msl libraries first then text libraries.
		/// </summary>
		/// <param name="sequence">Modified sequence in mzLib bracket notation.</param>
		/// <param name="charge">Precursor charge state.</param>
		/// <param name="librarySpectrum">
		///   Set to the matching spectrum on success; null on failure.
		/// </param>
		/// <returns>True if the spectrum was found; false otherwise.</returns>
		public bool TryGetSpectrum(string sequence, int charge, out LibrarySpectrum librarySpectrum)
		{
			// Check MSL libraries first — they maintain their own index and need no buffer
			foreach (var mslLib in _mslLibraries.Values)
			{
				if (mslLib.TryGetLibrarySpectrum(sequence, charge, out librarySpectrum))
					return true;
			}

			string lookupString = sequence + "/" + charge;
			librarySpectrum = null;

			// look up in buffer to see if this library spectrum was read in already
			if (LibrarySpectrumBuffer.TryGetValue(lookupString, out var spectrum))
			{
				librarySpectrum = spectrum;

				if (librarySpectrum.Name != lookupString)
				{
					throw new MzLibException("Bad spectral library formatting or indexing: Found \""
						+ librarySpectrum.Name + "\" but expected \"" + lookupString + "\"");
				}

				return true;
			}

			// go find the library spectrum in the spectral library file
			if (SequenceToFileAndLocation.TryGetValue(lookupString, out var value))
			{
				lock (StreamReaders[value.filePath])
				{
					librarySpectrum = ReadSpectrumFromLibraryFile(value.filePath, value.byteOffset);

					if (librarySpectrum.Name != lookupString)
					{
						throw new MzLibException("Bad spectral library formatting or indexing: Found \""
							+ librarySpectrum.Name + "\" but expected \"" + lookupString + "\"");
					}

					// add this item to the buffer
					lock (LibrarySpectrumBuffer)
					{
						lock (LibrarySpectrumBufferList)
						{
							LibrarySpectrumBuffer.TryAdd(lookupString, librarySpectrum);

							LibrarySpectrumBufferList.Enqueue(lookupString);


							// remove items from buffer if the buffer is at max capacity
							while (LibrarySpectrumBuffer.Count > MaxElementsInBuffer)
							{
								var item = LibrarySpectrumBufferList.Dequeue();
								LibrarySpectrumBuffer.Remove(item);
							}
						}
					}
				}

				return true;
			}

			return false;
		}

		public IEnumerable<LibrarySpectrum> GetAllLibrarySpectra()
		{
			foreach (var item in SequenceToFileAndLocation)
			{
				yield return ReadSpectrumFromLibraryFile(item.Value.filePath, item.Value.byteOffset);
			}
		}

		/// <summary>
		/// Closes all open <see cref="StreamReader"/> handles for text libraries and disposes
		/// all <see cref="MslLibrary"/> instances (which releases any open file handles held
		/// in index-only mode).  Safe to call multiple times.
		/// </summary>
		public void CloseConnections()
		{
			// Close text-library stream readers
			foreach (var item in StreamReaders)
			{
				item.Value.Close();
			}

			// Dispose MSL libraries — releases the FileStream in index-only mode
			foreach (var mslLib in _mslLibraries.Values)
			{
				mslLib.Dispose();
			}
		}

		private LibrarySpectrum ReadSpectrumFromLibraryFile(string path, long byteOffset)
		{
			if (!StreamReaders.TryGetValue(path, out var reader))
			{
				throw new MzLibException("????");
				//IndexSpectralLibrary(path);

				//if (!StreamReaders.TryGetValue(path, out reader))
				//{
				//    // TODO: throw an exception
				//}
			}

			// seek to the byte of the scan
			reader.BaseStream.Position = byteOffset;
			reader.DiscardBufferedData();

			// return the library spectrum
			if (path.Contains("pdeep"))
			{
				return ReadLibrarySpectrum_pDeep(reader);
			}
			else if (path.Contains("ms2pip"))
			{
				return ReadLibrarySpectrum_ms2pip(reader);
			}
			else
			{
				return ReadLibrarySpectrum(reader);
			}
		}

		private LibrarySpectrum ReadLibrarySpectrum(StreamReader reader, bool onlyReadHeader = false)
		{
			char[] nameSplit = new char[] { '/' };
			char[] mwSplit = new char[] { ':' };
			char[] commentSplit = new char[] { ' ', ':', '=' };
			char[] modSplit = new char[] { '=', '/' };
			char[] fragmentSplit = new char[] { '\t', '\"', ')', '/' };
			char[] neutralLossSplit = new char[] { '-' };

			bool readingPeaks = false;
			string sequence = null;
			int z = 2;
			double precursorMz = 0;
			double rt = 0;
			List<MatchedFragmentIon> matchedFragmentIons = new List<MatchedFragmentIon>();

			while (reader.Peek() > 0)
			{
				string line = reader.ReadLine();
				string[] split;

				if (line.StartsWith("Name", StringComparison.InvariantCultureIgnoreCase))
				{
					if (CrosslinkLibrarySpectrum.CrosslinkRegex.Match(line).Success)
					{
						return ReadLibrarySpectrum_Crosslink(reader, line, onlyReadHeader);
					}

					if (sequence != null)
					{
						return new LibrarySpectrum(sequence, precursorMz, z, matchedFragmentIons, rt);
					}

					split = line.Split(nameSplit);

					// get sequence
					sequence = split[0].Replace("Name:", string.Empty).Trim();

					// get charge
					z = int.Parse(split[1].Trim());
				}
				else if (line.StartsWith("MW", StringComparison.InvariantCultureIgnoreCase))
				{
					split = line.Split(mwSplit);

					// get precursor m/z
					precursorMz = double.Parse(split[1].Trim(), CultureInfo.InvariantCulture);
				}
				else if (line.StartsWith("Comment", StringComparison.InvariantCultureIgnoreCase))
				{
					split = line.Split(commentSplit);

					// get precursor m/z if not defined yet
					if (precursorMz == 0)
					{
						int indOfParent = Array.IndexOf(split, "Parent");
						if (indOfParent > 0)
						{
							precursorMz = double.Parse(split[indOfParent + 1], CultureInfo.InvariantCulture);
						}
					}

					// get RT
					int indOfRt = Array.IndexOf(split, "iRT");
					if (indOfRt > 0)
					{
						rt = double.Parse(split[indOfRt + 1], CultureInfo.InvariantCulture);
					}
					else
					{
						indOfRt = Array.IndexOf(split, "RT");

						if (indOfRt > 0)
						{
							rt = double.Parse(split[indOfRt + 1], CultureInfo.InvariantCulture);
						}
					}

					// get mods
					// be careful about spaces! mod names can have spaces in them
					StringBuilder sb = new StringBuilder();
					int ind = line.IndexOf("Mods", StringComparison.InvariantCultureIgnoreCase);

					if (ind > 0)
					{
						bool readingModName = false;
						int bracketCount = 0;

						for (int i = ind; i < line.Length; i++)
						{
							if (line[i] == ' ' && !readingModName)
							{
								break;
							}
							if (line[i] == '[')
							{
								bracketCount++;
								readingModName = true;
							}
							else if (line[i] == ']')
							{
								bracketCount--;

								if (bracketCount == 0)
								{
									readingModName = false;
								}
							}

							sb.Append(line[i]);
						}

						if (sb.ToString() != "Mods=0")
						{
							split = sb.ToString().Split(modSplit);

							for (int i = split.Length - 1; i >= 2; i--)
							{
								string modString = split[i];

								string[] modInfo = modString.Split(',');
								int modPosition = int.Parse(modInfo[0]);
								string modName = modInfo[2];
								string modNameNoBrackets = modName;

								if (modName.StartsWith('['))
								{
									modNameNoBrackets = modName.Substring(1, modName.Length - 2);
								}

								if (!ModificationConverter.AllKnownMods.Select(m => m.IdWithMotif).Contains(modNameNoBrackets))
								{
									if (PrositToMetaMorpheusModDictionary.TryGetValue(modName, out var metaMorpheusMod))
									{
										modName = metaMorpheusMod;
									}
								}

								// add the mod name into the sequence
								string leftSeq = sequence.Substring(0, modPosition + 1);
								string rightSeq = sequence.Substring(modPosition + 1);

								sequence = leftSeq + modName + rightSeq;
							}
						}
					}
				}
				else if (line.StartsWith("Num peaks", StringComparison.InvariantCultureIgnoreCase))
				{
					if (onlyReadHeader)
					{
						return new LibrarySpectrum(sequence, precursorMz, z, matchedFragmentIons, rt);
					}

					// this assumes that the peaks are listed after the "Num peaks" line
					readingPeaks = true;
				}
				else if (readingPeaks)
				{
					matchedFragmentIons.Add(ReadFragmentIon(line, fragmentSplit, neutralLossSplit, sequence));
				}
			}

			return new LibrarySpectrum(sequence, precursorMz, z, matchedFragmentIons, rt);
		}

		private LibrarySpectrum ReadLibrarySpectrum_pDeep(StreamReader reader, bool onlyReadHeader = false)
		{
			char[] nameSplit = new char[] { '/', '_' };
			char[] mwSplit = new char[] { ':' };
			char[] commentSplit = new char[] { ' ', ':', '=' };
			char[] modSplit = new char[] { '/', '(', ')' };
			char[] fragmentSplit = new char[] { '\t', '/' };

			bool readingPeaks = false;
			string sequence = null;
			int z = 2;
			double precursorMz = 0;
			double rt = 0;
			List<MatchedFragmentIon> matchedFragmentIons = new List<MatchedFragmentIon>();

			while (reader.Peek() > 0)
			{
				string line = reader.ReadLine();
				string[] split;

				if (line.StartsWith("Name", StringComparison.InvariantCultureIgnoreCase))
				{
					if (sequence != null)
					{
						return new LibrarySpectrum(sequence, precursorMz, z, matchedFragmentIons, rt);
					}

					split = line.Split(nameSplit);

					// get sequence
					sequence = split[0].Replace("Name:", string.Empty).Trim();

					// get charge
					z = int.Parse(split[1].Trim());

					string[] mods = split[2].Split(modSplit, StringSplitOptions.RemoveEmptyEntries);
					for (int i = mods.Length - 1; i > 0; i--)
					{
						string[] modInfo = mods[i].Split(',');
						int index = Convert.ToInt32(modInfo[0]);
						string mod = modInfo[2];
						string metaMorpheusMod = pDeepToMetaMorpheusModDictionary[mod];
						//add the mod into the sequence
						string leftSeq = sequence.Substring(0, index + 1);
						string rightSeq = sequence.Substring(index + 1);
						sequence = leftSeq + metaMorpheusMod + rightSeq;
					}

				}
				else if (line.StartsWith("Comment", StringComparison.InvariantCultureIgnoreCase))
				{
					split = line.Split(commentSplit);

					// get precursor m/z in comment
					int indOfParent = Array.IndexOf(split, "Parent");
					if (indOfParent > 0)
					{
						precursorMz = double.Parse(split[indOfParent + 1]);
					}

					// get RT
					int indOfRt = Array.IndexOf(split, "RTInSeconds");
					if (indOfRt > 0)
					{
						rt = double.Parse(split[indOfRt + 1]);
					}
				}
				else if (line.StartsWith("Num peaks", StringComparison.InvariantCultureIgnoreCase))
				{
					if (onlyReadHeader)
					{
						return new LibrarySpectrum(sequence, precursorMz, z, matchedFragmentIons, rt);
					}

					// this assumes that the peaks are listed after the "Num peaks" line
					readingPeaks = true;
				}
				else if (readingPeaks && line != "")
				{
					split = line.Split(fragmentSplit, StringSplitOptions.RemoveEmptyEntries);

					// read fragment m/z
					var experMz = double.Parse(split[0], CultureInfo.InvariantCulture);

					// read fragment intensity
					var experIntensity = double.Parse(split[1], CultureInfo.InvariantCulture);

					// read fragment type, number     
					string fragmentType = split[2].ToCharArray()[0].ToString();
					int fragmentNumber = int.Parse(new string(split[2].Split(new char[] { '^' })[0].Where(Char.IsDigit).ToArray()));
					int fragmentCharge = 1;

					if (split[2].Contains('^'))
					{
						fragmentCharge = int.Parse(split[2].Split('^')[1]);
					}
					ProductType peakProductType = (ProductType)Enum.Parse(typeof(ProductType), fragmentType, true);

					//TODO: figure out terminus
					FragmentationTerminus terminus = (FragmentationTerminus)Enum.Parse(typeof(FragmentationTerminus), "None", true);

					//TODO: figure out amino acid position
					var product = new Product(peakProductType, terminus, experMz, fragmentNumber, 0, 0);

					matchedFragmentIons.Add(new MatchedFragmentIon(product, experMz, experIntensity, fragmentCharge));
				}
			}

			return new LibrarySpectrum(sequence, precursorMz, z, matchedFragmentIons, rt);
		}

		private LibrarySpectrum ReadLibrarySpectrum_ms2pip(StreamReader reader, bool onlyReadHeader = false)
		{
			char[] nameSplit = new char[] { '/' };
			char[] mwSplit = new char[] { ':' };
			char[] commentSplit = new char[] { ' ', ':', '=' };
			char[] modSplit = new char[] { '=', '/' };
			char[] fragmentSplit = new char[] { '\t', '\"', ')', '/' };
			char[] neutralLossSplit = new char[] { '-' };

			bool readingPeaks = false;
			string sequence = null;
			int z = 2;
			double precursorMz = 0;
			double rt = 0;
			List<MatchedFragmentIon> matchedFragmentIons = new List<MatchedFragmentIon>();

			while (reader.Peek() > 0)
			{
				string line = reader.ReadLine();
				if (!line.IsNotNullOrEmpty())
				{
					continue;
				}
				string[] split;

				if (line.StartsWith("Name", StringComparison.InvariantCultureIgnoreCase))
				{
					if (CrosslinkLibrarySpectrum.CrosslinkRegex.Match(line).Success)
					{
						return ReadLibrarySpectrum_Crosslink(reader, line, onlyReadHeader);
					}

					if (sequence != null)
					{
						return new LibrarySpectrum(sequence, precursorMz, z, matchedFragmentIons, rt);
					}

					split = line.Split(nameSplit);

					// get sequence
					sequence = split[0].Replace("Name:", string.Empty).Trim();

					// get charge
					z = int.Parse(split[1].Trim());
				}
				else if (line.StartsWith("MW", StringComparison.InvariantCultureIgnoreCase))
				{
					continue;
				}
				else if (line.StartsWith("Comment", StringComparison.InvariantCultureIgnoreCase))
				{
					split = line.Split(commentSplit);

					// get precursor m/z if not defined yet
					if (precursorMz == 0)
					{
						int indOfParent = Array.IndexOf(split, "Parent");
						if (indOfParent > 0)
						{
							precursorMz = double.Parse(split[indOfParent + 1], CultureInfo.InvariantCulture);
						}
					}

					// get mods
					// be careful about spaces! mod names can have spaces in them
					StringBuilder sb = new StringBuilder();
					int ind = line.IndexOf("Mods", StringComparison.InvariantCultureIgnoreCase);

					if (ind > 0)
					{
						bool readingModName = false;
						int bracketCount = 0;

						for (int i = ind; i < line.Length; i++)
						{
							if (line[i] == ' ' && !readingModName)
							{
								break;
							}
							sb.Append(line[i]);
						}

						if (sb.ToString() != "Mods=0")
						{
							split = sb.ToString().Split(modSplit);

							for (int i = split.Length - 1; i >= 2; i--)
							{
								string modString = split[i];

								string[] modInfo = modString.Split(',');
								int modPosition = int.Parse(modInfo[0]);
								string modName = modInfo[2];
								string modNameNoBrackets = modName;

								if (modName.StartsWith('['))
								{
									modNameNoBrackets = modName.Substring(1, modName.Length - 2);
								}

								if (!ModificationConverter.AllKnownMods.Select(m => m.IdWithMotif).Contains(modNameNoBrackets))
								{
									if (PrositToMetaMorpheusModDictionary.TryGetValue(modName, out var metaMorpheusMod))
									{
										modName = metaMorpheusMod;
									}
								}

								// add the mod name into the sequence
								string leftSeq = sequence.Substring(0, modPosition);
								string rightSeq = sequence.Substring(modPosition);

								sequence = leftSeq + modName + rightSeq;
							}
						}
					}
				}
				else if (line.StartsWith("Num peaks", StringComparison.InvariantCultureIgnoreCase))
				{
					if (onlyReadHeader)
					{
						return new LibrarySpectrum(sequence, precursorMz, z, matchedFragmentIons, rt);
					}

					// this assumes that the peaks are listed after the "Num peaks" line
					readingPeaks = true;
				}
				else if (readingPeaks)
				{
					matchedFragmentIons.Add(ReadFragmentIon(line, fragmentSplit, neutralLossSplit, sequence));
				}
			}

			return new LibrarySpectrum(sequence, precursorMz, z, matchedFragmentIons, rt);
		}

		internal CrosslinkLibrarySpectrum ReadLibrarySpectrum_Crosslink(StreamReader reader, string nameLine, bool onlyReadHeader)
		{
			char[] nameSplit = new char[] { '/' };
			char[] mwSplit = new char[] { ':' };
			char[] commentSplit = new char[] { ' ', ':', '=' };
			char[] fragmentSplit = new char[] { '\t', '\"', ')', '/' };
			char[] neutralLossSplit = new char[] { '-' };

			string[] splitNameLine = nameLine.Split(nameSplit);
			string uniqueSequence = splitNameLine[0].Replace("Name:", string.Empty).Trim();
			string alphaSequence = "";
			string betaSequence = "";
			string[] splitAlphaBetaSequence = new Regex(pattern: @"\(\d+\)").Split(uniqueSequence);
			if (splitAlphaBetaSequence.Length >= 2)
			{
				alphaSequence = splitAlphaBetaSequence[0];
				betaSequence = splitAlphaBetaSequence[1];
			}
			else if (splitAlphaBetaSequence.Length == 1)
			{
				alphaSequence = splitAlphaBetaSequence[0];
			}
			int z = int.Parse(splitNameLine[1].Trim());
			double precursorMz = 0;
			double rt = 0;
			int indOfRt = -1;
			bool readingPeaks = false;
			List<MatchedFragmentIon> alphaPeptideIons = new List<MatchedFragmentIon>();
			List<MatchedFragmentIon> betaPeptideIons = new List<MatchedFragmentIon>();

			while (reader.Peek() > 0)
			{
				string line = reader.ReadLine();
				string[] split;

				if (line.StartsWith("Name", StringComparison.InvariantCultureIgnoreCase))
				{
					break;
				}
				else if (line.StartsWith("MW", StringComparison.InvariantCultureIgnoreCase))
				{
					split = line.Split(mwSplit);

					// get precursor m/z
					precursorMz = double.Parse(split[1].Trim(), CultureInfo.InvariantCulture);
				}
				else if (line.StartsWith("Comment", StringComparison.InvariantCultureIgnoreCase))
				{
					split = line.Split(commentSplit);

					// get precursor m/z if not defined yet
					if (precursorMz == 0)
					{
						int indOfParent = Array.IndexOf(split, "Parent");
						if (indOfParent > 0)
						{
							precursorMz = double.Parse(split[indOfParent + 1], CultureInfo.InvariantCulture);
						}
					}


					indOfRt = Array.IndexOf(split, "RT");

					if (indOfRt > 0)
					{
						rt = double.Parse(split[indOfRt + 1], CultureInfo.InvariantCulture);
					}
				}
				else if (line.StartsWith("Num alpha peaks", StringComparison.InvariantCultureIgnoreCase))
				{
					if (onlyReadHeader)
					{
						CrosslinkLibrarySpectrum betaPeptideSpectrumHeaderOnly = new CrosslinkLibrarySpectrum(uniqueSequence, precursorMz, z, betaPeptideIons, rt);
						return new CrosslinkLibrarySpectrum(uniqueSequence, precursorMz, z, alphaPeptideIons, rt, betaPeptideSpectrumHeaderOnly);
					}
					// this assumes that the peaks are listed after the "Num peaks" line
					readingPeaks = true;
				}
				else if (readingPeaks)
				{
					bool isBetaPeptideIon = line.Contains("BetaPeptide");
					string peptideSequence = isBetaPeptideIon ? betaSequence : alphaSequence;
					MatchedFragmentIon fragmentIon =
						ReadFragmentIon(line, fragmentSplit, neutralLossSplit, peptideSequence);
					if (isBetaPeptideIon)
					{
						betaPeptideIons.Add(fragmentIon);
					}
					else
					{
						alphaPeptideIons.Add(fragmentIon);
					}
				}
			}

			CrosslinkLibrarySpectrum betaPeptideSpectrum = new CrosslinkLibrarySpectrum(uniqueSequence, precursorMz, z, betaPeptideIons, rt);
			return new CrosslinkLibrarySpectrum(uniqueSequence, precursorMz, z, alphaPeptideIons, rt, betaPeptideSpectrum);
		}

		/// <summary>
		/// Parses a single peak line from an MSP-format spectral library into a
		/// <see cref="MatchedFragmentIon"/>.  Handles both standard terminal ions
		/// (e.g. <c>b5</c>, <c>y12^2</c>, <c>y8-17.0265</c>) and internal fragment ions
		/// (e.g. <c>bIb[3-6]</c>, <c>aIb[2-5]^2</c>) written by
		/// <see cref="Omics.SpectralMatch.MslSpectralLibrary.MslFragmentIon.Annotation"/>.
		/// Does not work with P-Deep libraries.
		/// </summary>
		/// <param name="fragmentIonLine">
		///   A single tab-delimited peak line containing m/z, intensity, and annotation fields.
		/// </param>
		/// <param name="fragmentSplit">
		///   Characters used to split the peak line into fields (typically tab, quote, paren, slash).
		/// </param>
		/// <param name="neutralLossSplit">
		///   Characters used to split the annotation field when extracting a neutral-loss mass
		///   (typically just the hyphen <c>-</c>).
		/// </param>
		/// <param name="peptideSequence">
		///   The unmodified amino-acid sequence of the parent peptide, used to compute
		///   <c>ResiduePosition</c> for terminus-specific ions.  May be null or empty, in which
		///   case a default peptide length of 25 is assumed.
		/// </param>
		/// <returns>A fully populated <see cref="MatchedFragmentIon"/>.</returns>
		public static MatchedFragmentIon ReadFragmentIon(string fragmentIonLine, char[] fragmentSplit,
			char[] neutralLossSplit, string peptideSequence)
		{
			// Split the line into columns: [0] = m/z, [1] = intensity, [2] = annotation
			string[] split = fragmentIonLine.Split(fragmentSplit, StringSplitOptions.RemoveEmptyEntries);

			// Column 0: observed fragment m/z
			var experMz = double.Parse(split[0], CultureInfo.InvariantCulture);

			// Column 1: relative or absolute fragment intensity
			var experIntensity = double.Parse(split[1], CultureInfo.InvariantCulture);

			// Column 2: annotation string — try internal ion regex first, then terminal regex
			string annotation = split[2];

			// ── Internal ion path ─────────────────────────────────────────────────────
			// Internal ions are annotated as {PrimaryType}I{SecondaryType}[{start}-{end}]
			// e.g. "bIb[3-6]" or "aIb[2-5]^2". The standard IonParserRegex cannot match
			// the "I" separator or the bracket-enclosed residue range.
			Match internalMatch = InternalIonRegex.Match(annotation);
			if (internalMatch.Success)
			{
				// Group 1: primary product type string, e.g. "b" or "a"
				string primaryTypeStr = internalMatch.Groups[1].Value;

				// Group 2: secondary product type string (C-terminal boundary), e.g. "b" or "y"
				string secondaryTypeStr = internalMatch.Groups[2].Value;

				// Group 3: start residue number (N-terminal boundary of the internal span)
				int startResidue = int.Parse(internalMatch.Groups[3].Value);

				// Group 4: end residue number (C-terminal boundary of the internal span)
				int endResidue = int.Parse(internalMatch.Groups[4].Value);

				// Group 5: optional charge; default to 1 when absent
				int internalCharge = internalMatch.Groups[5].Success
					? int.Parse(internalMatch.Groups[5].Value)
					: 1;

				// Parse product types — throws ArgumentException for unknown type names
				ProductType primaryType = (ProductType)Enum.Parse(typeof(ProductType), primaryTypeStr, true);
				ProductType secondaryType = (ProductType)Enum.Parse(typeof(ProductType), secondaryTypeStr, true);

				// Internal ions have no peptide terminus association.
				// fragmentNumber holds the N-terminal boundary (startResidue) and
				// secondaryFragmentNumber holds the C-terminal boundary (endResidue),
				// consistent with how MslFragmentIon and FromLibrarySpectrum represent
				// internal ions in the MSL binary format.
				var internalProduct = new Product(
					primaryType,
					FragmentationTerminus.None,
					experMz.ToMass(internalCharge),
					fragmentNumber: startResidue,
					residuePosition: startResidue,
					neutralLoss: 0,
					secondaryProductType: secondaryType,
					secondaryFragmentNumber: endResidue);

				return new MatchedFragmentIon(internalProduct, experMz, experIntensity, internalCharge);
			}

			// ── Terminal ion path (unchanged from original) ───────────────────────────
			Match regexMatchResult = IonParserRegex.Match(annotation);

			// Group 1: ion type letters (e.g. "b", "y", "zDot")
			string fragmentType = regexMatchResult.Groups[1].Value;

			// Group 2: fragment series number
			int fragmentNumber = int.Parse(regexMatchResult.Groups[2].Value);

			// Group 3: charge digit(s) after "^", or empty
			int fragmentCharge = 1;
			if (regexMatchResult.Groups.Count > 3 && !string.IsNullOrWhiteSpace(regexMatchResult.Groups[3].Value))
			{
				fragmentCharge = int.Parse(regexMatchResult.Groups[3].Value);
			}

			// Neutral-loss mass: extracted from a trailing "-{mass}" in the annotation
			double neutralLoss = 0;
			if (annotation.Contains("-") && fragmentCharge > 0)
			{
				string[] neutralLossInformation = annotation.Split(neutralLossSplit, StringSplitOptions.RemoveEmptyEntries).ToArray();
				neutralLoss = double.Parse(neutralLossInformation[1]);
			}
			if (fragmentCharge < 0)
			{
				string[] neutralLossInformation = annotation.Split(neutralLossSplit, StringSplitOptions.RemoveEmptyEntries).ToArray();
				if (neutralLossInformation.Length > 2)
					neutralLoss = double.Parse(neutralLossInformation[2]);
			}

			ProductType peakProductType = (ProductType)Enum.Parse(typeof(ProductType), fragmentType, true);

			// Default product for ProductTypes not in ProductTypeToFragmentationTerminus (e.g. "M" ions)
			Product product = new Product(peakProductType, (FragmentationTerminus)Enum.Parse(typeof(FragmentationTerminus),
				"None", true), experMz, fragmentNumber, 0, 0);

			if (TerminusSpecificProductTypes.ProductTypeToFragmentationTerminus.TryGetValue(peakProductType,
					out var terminus))
			{
				// peptideLength is used to compute ResiduePosition for C-terminal ions
				int peptideLength = peptideSequence.IsNotNullOrEmptyOrWhiteSpace() ? peptideSequence.Length : 25; // Arbitrary default peptide length
				product = new Product(peakProductType, terminus, experMz.ToMass(fragmentCharge), fragmentNumber,
					residuePosition: terminus == FragmentationTerminus.N ? fragmentNumber : peptideLength - fragmentNumber,
					neutralLoss);
			}

			return new MatchedFragmentIon(product, experMz, experIntensity, fragmentCharge);
		}

		private void IndexSpectralLibrary(string path)
		{
			var reader = new StreamReader(path);
			StreamReaders.Add(path, reader);

			reader.BaseStream.Position = 0;
			reader.DiscardBufferedData();

			while (reader.Peek() > 0)
			{
				long byteOffset = TextFileReading.GetByteOffsetAtCurrentPosition(reader);
				var line = reader.ReadLine().Trim();

				if (line.StartsWith("name", StringComparison.InvariantCultureIgnoreCase))
				{
					// seek back to beginning of line so the parser can read the "name" line
					reader.BaseStream.Position = byteOffset;
					reader.DiscardBufferedData();

					// parse the header
					//TODO: This is not very reliable way of determining the input library format but it works for our purpose and has no problem reading in libraries generated from MetaMorpheus.
					// A user-defined library format might be necessary if we want to use outputs from various external tools.
					LibrarySpectrum libraryItem;
					if (path.Contains("pdeep"))
					{
						libraryItem = ReadLibrarySpectrum_pDeep(reader, onlyReadHeader: true);
					}
					else if (path.Contains("ms2pip"))
					{
						libraryItem = ReadLibrarySpectrum_ms2pip(reader, onlyReadHeader: true);
					}
					else
					{
						libraryItem = ReadLibrarySpectrum(reader, onlyReadHeader: true);
					}

					// add the spectrum to the index
					SequenceToFileAndLocation.TryAdd(libraryItem.Name, (path, byteOffset));
				}
			}
		}
	}
}