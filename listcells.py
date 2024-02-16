
L1_NGC = [f"L1_NGC-DA_bNAC219_{i}" for i in range(1,6)]
L23 = [f"L23_PC_cADpyr229_{i}" for i in range(1,6)]
L4_PC = [f"L4_PC_cADpyr230_{i}" for i in range(1,6)]
L4_LBC = [f"L4_LBC_cACint209_{i}" for i in range(1,6)]
L5 = [f"L5_TTPC2_cADpyr232_{i}" for i in range(1,6)]
L6_TPC = [f"L6_TPC_L4_cADpyr231_{i}" for i in range(1,6)]

L6 = ["L6_TPC_L1_cADpyr231_1","L6_TPC_L4_cADpyr231_1", "L6_UTPC_cADpyr231_1"]

SST = ['L4_MC_cACint209_1','L4_MC_cACint209_2','L4_MC_cACint209_3','L4_MC_cACint209_4','L4_MC_cACint209_5']

VIP = ['L4_BP_bNAC219_1','L4_BP_bNAC219_2', 'L4_BP_bNAC219_3', 'L4_BP_bNAC219_4', 'L4_BP_bNAC219_5',
       'L4_BP_cACint209_1', 'L4_BP_cACint209_2', 'L4_BP_cACint209_3','L4_BP_cACint209_4', 'L4_BP_cACint209_5' ]

PV = ['L4_LBC_dNAC222_1', 'L4_NBC_cNAC187_1', 'L4_NBC_dNAC222_1', 
      'L4_SBC_bNAC219_1',  'L4_SBC_cACint209_1','L4_SBC_dNAC222_1',
      'L4_ChC_cACint209_1','L4_ChC_cNAC187_1', 'L4_ChC_dNAC222_1']
PC = L23+L5+L6_TPC

allCells = L5  + L23+ L4_LBC  + L1_NGC + L6_TPC + VIP + SST+PV

L4_inhib = [  
'L4_BP_bNAC219_1','L4_BP_cACint209_1',
'L4_ChC_cACint209_1','L4_ChC_cNAC187_1', 'L4_ChC_dNAC222_1',
'L4_LBC_dNAC222_1', 
'L4_MC_cACint209_1',
'L4_NBC_cNAC187_1', 'L4_NBC_dNAC222_1', 
'L4_SBC_bNAC219_1',  'L4_SBC_cACint209_1','L4_SBC_dNAC222_1']


alpha10Hz = {
	"L1_NGC-DA_bNAC219_1" : 17.5,
    "L1_NGC-DA_bNAC219_2" : 45,
    "L1_NGC-DA_bNAC219_3" : 30,
    "L1_NGC-DA_bNAC219_4" : 13,
    "L1_NGC-DA_bNAC219_5" : 3.2,

	"L4_BP_cACint209_1" : 135 ,
	"L4_BP_cACint209_2" : 45 ,
	"L4_BP_cACint209_3" : 130 ,
	"L4_BP_cACint209_4" : 65 ,
	"L4_BP_cACint209_5" : 100 ,
 
	"L4_BP_bNAC219_1" : 4   , 
 	"L4_BP_bNAC219_2" : 12	,
	"L4_BP_bNAC219_3" : 12.5,
	"L4_BP_bNAC219_4" : 4.6,
	"L4_BP_bNAC219_5" : 3.7,
 
	"L4_LBC_dNAC222_1" : 0.4315	,
    "L4_LBC_dNAC222_2" : 0.31	,
    "L4_LBC_dNAC222_3" : 0.44	,
    "L4_LBC_dNAC222_4" : 0.42	,
    "L4_LBC_dNAC222_5" : 0.43	,    
	"L4_LBC_cACint209_1":0.48,
	"L4_LBC_cACint209_2":0.70,
	"L4_LBC_cACint209_3":0.42,
	"L4_LBC_cACint209_4":0.58,
	"L4_LBC_cACint209_5":0.64,
 
	"L4_SBC_cACint209_1" : 0.62 ,
	"L4_SBC_cACint209_2" : 0.46 ,
	"L4_SBC_cACint209_3" : 0.78 ,
	"L4_SBC_cACint209_4" : 0.69 ,
	"L4_SBC_cACint209_5" : 0.47 ,
    "L4_SBC_bNAC219_1" :1.25,
	"L4_SBC_dNAC222_1" :0.9,
 
	"L4_NBC_cNAC187_1" : 0.47  ,
	"L4_NBC_dNAC222_1" : 0.27,
 
	"L4_ChC_dNAC222_1" : 0.76,
    "L4_ChC_cNAC187_1" : 0.474	,
	"L4_ChC_cACint209_1" : 0.405,
	"L4_ChC_cACint209_2" : 0.54,
	"L4_ChC_cACint209_3" : 0.45,
	"L4_ChC_cACint209_4" : 0.42,
	"L4_ChC_cACint209_5" : 0.75,
    	
	"L4_MC_cACint209_1" : 1.10,
	"L4_MC_cACint209_2" : 3.25,
	"L4_MC_cACint209_3" : 0.7,
 	"L4_MC_cACint209_4" : 0.95,
  	"L4_MC_cACint209_5" : 1.1,
 
	"L23_PC_cADpyr229_1" : 4.5	,
	"L23_PC_cADpyr229_2": 4.8,
	"L23_PC_cADpyr229_3" : 6.1,	
	"L23_PC_cADpyr229_4" : 5.8,	
	"L23_PC_cADpyr229_5" : 3,
    
	"L4_PC_cADpyr230_1" : 0.78,
	"L4_PC_cADpyr230_2" : 3.2,
	"L4_PC_cADpyr230_3" : 2.1,
	"L4_PC_cADpyr230_4" : 3.1,
	"L4_PC_cADpyr230_5" : 2.5,
    
	"L5_TTPC2_cADpyr232_1" : 0.61,
	"L5_TTPC2_cADpyr232_2" : 0.642,	
	"L5_TTPC2_cADpyr232_3" : 0.617	,
	"L5_TTPC2_cADpyr232_4" : 1.23	,
	"L5_TTPC2_cADpyr232_5" : 0.772	,
 
	"L5_TTPC1_cADpyr232_1" : 0.615,
	"L5_UTPC_cADpyr232_1": 1.15,
    
    "L6_TPC_L4_cADpyr231_1": 1.5,
    "L6_TPC_L4_cADpyr231_2": 1.5,
    "L6_TPC_L4_cADpyr231_3": 1.9,
    "L6_TPC_L4_cADpyr231_4": 1.3,
    "L6_TPC_L4_cADpyr231_5": 1.7,
    
	"L6_TPC_L1_cADpyr231_1": 1.45,
    "L6_UTPC_cADpyr231_1": 1.7,
}