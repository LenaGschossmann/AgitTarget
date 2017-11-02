#include <stdlib.h>
#include <math.h>
#ifdef MEX
#include "mex.h"
#else
#define mxCalloc calloc
#define mxFree free
#endif
#include "result.h"
#include "trials.h"

#ifdef MEX
#include "mex.h"
char * mxstr(const mxArray *string_array_ptr)
/* convert Matlab string to C char pointer */
{
    char *buf;
    int buflen;

    buflen = mxGetNumberOfElements(string_array_ptr) + 1;
    buf = (char *)mxCalloc(buflen, sizeof(char));
    if (mxGetString(string_array_ptr, buf, buflen) != 0)
        return NULL;
    return buf;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* i/o validation */
    if (nrhs < 1)
        mexErrMsgTxt("Invalid Input: requires filename (string)");
    char * filename = mxstr(prhs[0]);
    if (!filename)
        mexErrMsgTxt("Input argument is not a string.");
    char ** filter = (char**)mxCalloc(nrhs-1, sizeof(char*));
    for (int i=1; i<nrhs; i++)
        filter[i-1] = mxstr(prhs[i]);

    /* first pass: get trial counts*/
    Trials * trials = new Trials(filter, nrhs-1);
    if (! trials->readfile(filename))
        mexErrMsgTxt("Cannot read input file");

    plhs[0] = trials->data_struct();
    if (nlhs == 2)
        plhs[1] = trials->info_struct();
    delete trials;
    mxFree(filter);
    if (!plhs[0])
        mexErrMsgTxt("Cannot read input file on second pass.");

    return;
}

#else
int main(int argc, char * argv[])
{
    int filtermsgcount = 2;
    char * filtermsg[2] = { "APLAY", "FILL" };
    cout << "DEBUG" << endl;
    Trials t(filtermsg,filtermsgcount);
    t.readfile(argv[1]);
}
#endif
