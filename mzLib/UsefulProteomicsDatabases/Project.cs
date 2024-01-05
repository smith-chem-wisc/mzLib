using System;
using System.Net.Http;
using System.Threading.Tasks;
using Newtonsoft.Json;

public class Project
{
    private static readonly HttpClient client = new HttpClient();
    private const string ApiBaseUrl = "https://www.ebi.ac.uk/pride/ws/archive/v2/";
    private const string PrivateApiBaseUrl = "https://www.ebi.ac.uk/pride/private/ws/archive/v2/";

    public async Task<string> GetProjectsAsync(int pageSize, int page, string sortDirection, string sortConditions)
    {
        var requestUrl = $"{ApiBaseUrl}projects?pageSize={pageSize}&page={page}&sortDirection={sortDirection}&sortConditions={sortConditions}";
        var response = await client.GetStringAsync(requestUrl);
        return response;
    }

    public static async Task<string> GetReanalysisProjectsByAccessionAsync(string accession)
    {
        var requestUrl = $"{ApiBaseUrl}projects/reanalysis/{accession}";
        //var headers = "Accept": "application/JSON";
        var response = await client.GetStringAsync(requestUrl);
        return response;
    }

    public async Task<string> GetByAccessionAsync(string accession)
    {
        var requestUrl = $"{ApiBaseUrl}projects/{accession}";
        var response = await client.GetStringAsync(requestUrl);
        return response;
    }

    // Additional methods (GetFilesByAccession, GetPrivateFilesByAccession, etc.) would follow a similar pattern.
}