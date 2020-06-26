function xcp = CasADi_test_Speedyxcp

xcp.events     =  repmat(struct('id',{}, 'sampletime', {}, 'offset', {}), getNumEvents, 1 );
xcp.parameters =  repmat(struct('symbol',{}, 'size', {}, 'dtname', {}, 'baseaddr', {}), getNumParameters, 1 );
xcp.signals    =  repmat(struct('symbol',{}), getNumSignals, 1 );
xcp.models     =  cell(1,getNumModels);
    
xcp.models{1} = 'CasADi_test_Speedy';
   
   
         
xcp.events(1).id         = 0;
xcp.events(1).sampletime = 0.002;
xcp.events(1).offset     = 0.0;
    
xcp.signals(1).symbol =  'CasADi_test_Speedy_B.Constant';
    
xcp.signals(2).symbol =  'CasADi_test_Speedy_B.SFunction';
         
xcp.parameters(1).symbol = 'CasADi_test_Speedy_P.Constant_Value';
xcp.parameters(1).size   =  16;       
xcp.parameters(1).dtname = 'real_T'; 
xcp.parameters(2).baseaddr = '&CasADi_test_Speedy_P.Constant_Value[0]';     
         
xcp.parameters(2).symbol = 'CasADi_test_Speedy_P.SFunction_P1';
xcp.parameters(2).size   =  8;       
xcp.parameters(2).dtname = 'real_T'; 
xcp.parameters(3).baseaddr = '&CasADi_test_Speedy_P.SFunction_P1[0]';     
         
xcp.parameters(3).symbol = 'CasADi_test_Speedy_P.SFunction_P2';
xcp.parameters(3).size   =  1;       
xcp.parameters(3).dtname = 'real_T'; 
xcp.parameters(4).baseaddr = '&CasADi_test_Speedy_P.SFunction_P2';     

function n = getNumParameters
n = 3;

function n = getNumSignals
n = 2;

function n = getNumEvents
n = 1;

function n = getNumModels
n = 1;

