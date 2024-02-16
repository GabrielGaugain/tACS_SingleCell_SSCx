#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _Ca_reg();
extern void _Ca_HVA_reg();
extern void _Ca_LVAst_reg();
extern void _CaDynamics_E2_reg();
extern void _Ih_reg();
extern void _Im_reg();
extern void _K_Pst_reg();
extern void _K_Tst_reg();
extern void _KdShu2007_reg();
extern void _Nap_Et2_reg();
extern void _NaTa_t_reg();
extern void _NaTs2_t_reg();
extern void _netstimm_reg();
extern void _ProbAMPANMDA_EMS_reg();
extern void _ProbGABAAB_EMS_reg();
extern void _SK_E2_reg();
extern void _SKv3_1_reg();
extern void _vecevent_reg();
extern void _xtra_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," Ca.mod");
fprintf(stderr," Ca_HVA.mod");
fprintf(stderr," Ca_LVAst.mod");
fprintf(stderr," CaDynamics_E2.mod");
fprintf(stderr," Ih.mod");
fprintf(stderr," Im.mod");
fprintf(stderr," K_Pst.mod");
fprintf(stderr," K_Tst.mod");
fprintf(stderr," KdShu2007.mod");
fprintf(stderr," Nap_Et2.mod");
fprintf(stderr," NaTa_t.mod");
fprintf(stderr," NaTs2_t.mod");
fprintf(stderr," netstimm.mod");
fprintf(stderr," ProbAMPANMDA_EMS.mod");
fprintf(stderr," ProbGABAAB_EMS.mod");
fprintf(stderr," SK_E2.mod");
fprintf(stderr," SKv3_1.mod");
fprintf(stderr," vecevent.mod");
fprintf(stderr," xtra.mod");
fprintf(stderr, "\n");
    }
_Ca_reg();
_Ca_HVA_reg();
_Ca_LVAst_reg();
_CaDynamics_E2_reg();
_Ih_reg();
_Im_reg();
_K_Pst_reg();
_K_Tst_reg();
_KdShu2007_reg();
_Nap_Et2_reg();
_NaTa_t_reg();
_NaTs2_t_reg();
_netstimm_reg();
_ProbAMPANMDA_EMS_reg();
_ProbGABAAB_EMS_reg();
_SK_E2_reg();
_SKv3_1_reg();
_vecevent_reg();
_xtra_reg();
}
