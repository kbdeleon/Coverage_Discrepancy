/*
A KBase module: Coverage_Discrepancy
*/

module Coverage_Discrepancy {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */
    funcdef run_Coverage_Discrepancy(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;

};
