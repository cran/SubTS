#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .C calls */
extern void rDickman(void *, void *, void *);
extern void rF1(void *, void *, void *, void *);
extern void rF2(void *, void *, void *, void *);
extern void rPGamma(void *, void *, void *, void *, void *, void *, void *);
extern void rPRDTS(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void rPRDTSneg(void *, void *, void *, void *, void *);
extern void rTrunGamma(void *, void *, void *, void *, void *);
extern void rTrunSOptim(void *, void *, void *, void *, void *, void *, void *);
extern void rTrunTS(void *, void *, void *, void *, void *, void *, void *, void *);
extern void simCondS(void *, void *, void *);
extern void simTandW(void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"rDickman",    (DL_FUNC) &rDickman,     3},
    {"rF1",         (DL_FUNC) &rF1,          4},
    {"rF2",         (DL_FUNC) &rF2,          4},
    {"rPGamma",     (DL_FUNC) &rPGamma,      7},
    {"rPRDTS",      (DL_FUNC) &rPRDTS,      10},
    {"rPRDTSneg",   (DL_FUNC) &rPRDTSneg,    5},
    {"rTrunGamma",  (DL_FUNC) &rTrunGamma,   5},
    {"rTrunSOptim", (DL_FUNC) &rTrunSOptim,  7},
    {"rTrunTS",     (DL_FUNC) &rTrunTS,      8},
    {"simCondS",    (DL_FUNC) &simCondS,     3},
    {"simTandW",    (DL_FUNC) &simTandW,     4},
    {NULL, NULL, 0}
};

void R_init_SubTS(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
