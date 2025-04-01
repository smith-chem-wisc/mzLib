using Chemistry;
using CsvHelper;
using Omics.Fragmentation;
using Omics.SpectrumMatch;

namespace Readers.SpectrumLibraries
{
    public class MsFraggerSpeclibFile : SpectrumLibraryFile
    {
        public override SupportedFileType FileType { get; }
        public override Software Software { get; set; }

        public List<MsFraggerSpeclib> OriginalRecords { get; private set; } //single line in the .speclib file. combine ultiple lines to create a single library spectrum


        public MsFraggerSpeclibFile() : base()
        {
        }

        public MsFraggerSpeclibFile(string filePath) : base(filePath, Software.MsFragger)
        {
        }

        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), MsFraggerSpeclib.CsvConfiguration);
            OriginalRecords = csv.GetRecords<MsFraggerSpeclib>().ToList();
            List<LibrarySpectrum> librarySpectra = new List<LibrarySpectrum>();
            //each record contains an individual fragment ion. combine multiple records to create a single library spectrum
            //group by modified peptide and precursor charge state
            foreach (var spectrumFragmentGroup in OriginalRecords.GroupBy(p => new { p.ModifiedPeptide, p.PrecursorCharge }))
            {
                string sequence = spectrumFragmentGroup.Key.ModifiedPeptide;
                double precursorMz = spectrumFragmentGroup.First().PrecursorMz;
                int chargeState = spectrumFragmentGroup.Key.PrecursorCharge;
                List<MatchedFragmentIon> matchedFragmentIons = new();
                foreach (var spectrumFragment in spectrumFragmentGroup)
                {
                    double fragmentMz = spectrumFragment.ProductMz;
                    double intensity = spectrumFragment.LibraryIntensity;
                    int charge = spectrumFragment.FragmentCharge;

                    double neutralMass = fragmentMz.ToMass(charge);
                    ProductType productType = new ProductType();
                    switch (spectrumFragment.FragmentType)
                    {
                        case "a":
                            productType = ProductType.a;
                            break;
                        //case "a*":
                        //    pt = ProductType.aStar;
                        //    break;
                        //case "a°":
                        //    pt = ProductType.aDegree;
                        //    break;
                        //case "a-H2O":
                        //    pt = ProductType.aWaterLoss;
                        //    break;
                        //case "a-NH3":
                        //    pt = ProductType.aAmmoniaLoss;
                        //    break;
                        case "b":
                            productType = ProductType.b;
                            break;
                        //case "b-H2O":
                        //    pt = ProductType.bWaterLoss;
                        //    break;
                        //case "b-NH3":
                        //    pt = ProductType.bAmmoniaLoss;
                            break;
                        case "c":
                            productType = ProductType.c;
                            break;
                        //case "c-H2O":
                        //    pt = ProductType.cWaterLoss;
                        //    break;
                        //case "c-NH3":
                        //    pt = ProductType.cAmmoniaLoss;
                        //    break;
                        //case "d":
                        //    pt = ProductType.d;
                        //    break;
                        //case "d-H2O":
                        //    pt = ProductType.dWaterLoss;
                        //    break;
                        //case "d-NH3":
                        //    pt = ProductType.dAmmoniaLoss;
                        //    break;
                        //case "w":
                        //    pt = ProductType.w;
                        //    break;
                        //case "w-H2O":
                        //    pt = ProductType.wWaterLoss;
                        //    break;
                        //case "w-NH3":
                        //    pt = ProductType.wAmmoniaLoss;
                        //    break;
                        case "x":
                            productType = ProductType.x;
                            break;
                        //case "x-H2O":
                        //    pt = ProductType.xWaterLoss;
                        //    break;
                        //case "x-NH3":
                        //    pt = ProductType.xAmmoniaLoss;
                        //    break;
                        case "y":
                            productType = ProductType.y;
                            break;
                        //case "y-H2O":
                        //    pt = ProductType.yWaterLoss;
                        //    break;
                        //case "y-NH3":
                        //    pt = ProductType.yAmmoniaLoss;
                        //    break;
                        case "z":
                            productType = ProductType.z;
                            break;
                        //case "z-H2O":
                        //    pt = ProductType.zWaterLoss;
                        //    break;
                        //case "z-NH3":
                        //    pt = ProductType.zAmmoniaLoss;
                    }


                    double neutralLoss = 0;
                    //TODO: add corrct neutral losses
                    switch (spectrumFragment.FragmentLossType)
                    {
                        case "H2O":
                            neutralLoss = 18;
                            break;
                        case "NH3":
                            neutralLoss = 17;
                            break;
                    }

                    FragmentationTerminus fragmentationTerminus = new();

                    switch(productType)
                    {
                        case ProductType.a:
                        case ProductType.b:
                        case ProductType.c:
                            fragmentationTerminus = FragmentationTerminus.N;
                            break;
                        case ProductType.x:
                        case ProductType.y:
                        case ProductType.z:
                            fragmentationTerminus = FragmentationTerminus.C;
                            break;
                    }

                    int fragmentNumber = spectrumFragment.FragmentSeriesNumber;
                    int residuePosition = spectrumFragment.FragmentSeriesNumber;
                    Product neutralTheoreticalProduct = new Product(productType, fragmentationTerminus, neutralMass, fragmentNumber, residuePosition, neutralLoss);

                    MatchedFragmentIon matchedFragmentIon = new MatchedFragmentIon(neutralTheoreticalProduct, fragmentMz, intensity, charge);
                    matchedFragmentIons.Add(matchedFragmentIon);
                }
                double retentionTime = spectrumFragmentGroup.First().Tr_recalibrated;
                bool isDecoy = false;
                if(spectrumFragmentGroup.First().decoy.ToString() == "1")
                {
                    isDecoy = true;
                }

                librarySpectra.Add(new LibrarySpectrum(sequence, precursorMz, chargeState, matchedFragmentIons, retentionTime, isDecoy));
            }
            Results = librarySpectra;
        }

        public override void WriteResults(string outputPath)
        {
            if (!CanRead(outputPath))
                outputPath += FileType.GetFileExtension();

            using (var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), MsFraggerSpeclib.CsvConfiguration))
            {
                csv.WriteHeader<MsFraggerSpeclib>();
                foreach (var result in OriginalRecords)
                {
                    csv.NextRecord();
                    csv.WriteRecord(result);
                }
            }
        }
    }

}
