using System;
using System.Collections.Generic;
using System.IO;
using System.Text;

    public class FileCategory
    {
       public string Type { get; set; }
        public string CvLabel { get; set; }
        public string Accession { get; set; }
        public string Name { get; set; }
        public string Value { get; set; }
    }

    public class PublicFileLocation
    {
        public string Type { get; set; }
        public string CvLabel { get; set; }
        public string Accession { get; set; }
        public string Name { get; set; }
        public string Value { get; set; }
    }

    public class Link
    {
        public string Rel { get; set; }
        public string Href { get; set; }
    
    }

    public class _Link
    {
   
        public Link Self { get; set; }
        public Link Next { get; set; }
        public Link Previous { get; set; }
        public Link First { get; set; }
        public Link Last { get; set; }
        public Link Facets { get; set;}
        public Link DatasetFtpUrl { get; set;}
    }


    public class DetailType
    {
        public string type { get; set; }
        public string cvlabel { get; set; }
        public string accession { get; set; }
        public string name { get; set; }
    }

    public class Highlights
    {
        public string SampleAttributes { get; set; }
        public string Organisms { get; set; }
        public string ProjectDescription { get; set; }
        public string DataProcessingProtocol { get; set; }
        public string Organisms_facet { get; set; }
    }


    public class PRIDEEntry
    {
        public List<string> ProjectAccessions { get; set; }
        public string Accession { get; set; }
        public FileCategory FileCategory { get; set; }
        public string Checksum { get; set; }
        public List<PublicFileLocation> PublicFileLocations { get; set; }
        public long FileSizeBytes { get; set; }
        public string FileName { get; set; }
        public bool Compress { get; set; }
        public long SubmissionDate { get; set; }
        public long PublicationDate { get; set; }
        public long UpdatedDate { get; set; }
        public List<object> AdditionalAttributes { get; set; }
        public List<Link> Links { get; set; }
        public string title { get; set; }
        public string projectDescription { get; set; }
        public string sampleProcessingProtocol { get; set; }
        public string dataProcessingProtocol { get; set; }
        public List<string> KeyWords { get; set; }
    }

    public class PrideProject
    {
        public Highlights Highlights { get; set; }
        public string Accession { get; set; }
        public string title { get; set; }
        public List<object> AdditionalAttributes { get; set; }
        public string projectDescription { get; set; }
        public string sampleProcessingProtocol { get; set; }
        public string dataProcessingProtocol { get; set; }
        public List<string> KeyWords { get; set; }
        public string SubmissionDate { get; set; }
        public string PublicationDate { get; set; }
        //public List<string> submitters { get; set; }
        //public List<string> labPIs { get; set; }
        public List<string> affiliations { get; set; }
        public FileCategory instrument { get; set; }
        public FileCategory ExperimentType { get; set; }
       // public DetailType QuantificationMethods { get; set; }

        public List<FileCategory> Organisms { get; set; }
        //public String Organisms { get; set; }
        public List<FileCategory> OrganismParts { get; set; }
        //public String OrganismParts { get; set; }
        public List<FileCategory> Diseases { get; set; }
        public List<FileCategory> identifiedPTMStrings  { get; set; }
        public int QueryScore { get; set; }
        public _Link _Link { get; set; }
    }

    public class Page
    {
        public int Size;
        public int TotalElememts;
        public int TotalPages;
        public int Number;
    }

    public class _embedded
    {
        public List<PrideProject> Compactprojects { get; set; }
    }

    public class AllProjects_Embedded
    {
        public _embedded _Embedded { get; set; }
        public Page Page { get; set; }
        public _Link _Links { get; set; }
    }

   

