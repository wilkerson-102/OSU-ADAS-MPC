function pt=CasADi_test_Speedypt
pt = [];

  
pt(1).blockname = 'Constant';
pt(1).paramname = 'Value';
pt(1).class     = 'vector';
pt(1).nrows     = 16;
pt(1).ncols     = 1;
pt(1).subsource = 'SS_DOUBLE';
pt(1).ndims     = '2';
pt(1).size      = '[]';
pt(1).isStruct  = false;
pt(1).symbol     = 'CasADi_test_Speedy_P.Constant_Value';
pt(1).baseaddr   = '&CasADi_test_Speedy_P.Constant_Value[0]';
pt(1).dtname     = 'real_T';

pt(getlenPT) = pt(1);


  
pt(2).blockname = 'S-Function';
pt(2).paramname = 'P1';
pt(2).class     = 'vector';
pt(2).nrows     = 1;
pt(2).ncols     = 8;
pt(2).subsource = 'SS_DOUBLE';
pt(2).ndims     = '2';
pt(2).size      = '[]';
pt(2).isStruct  = false;
pt(2).symbol     = 'CasADi_test_Speedy_P.SFunction_P1';
pt(2).baseaddr   = '&CasADi_test_Speedy_P.SFunction_P1[0]';
pt(2).dtname     = 'real_T';



  
pt(3).blockname = 'S-Function';
pt(3).paramname = 'P2';
pt(3).class     = 'scalar';
pt(3).nrows     = 1;
pt(3).ncols     = 1;
pt(3).subsource = 'SS_DOUBLE';
pt(3).ndims     = '2';
pt(3).size      = '[]';
pt(3).isStruct  = false;
pt(3).symbol     = 'CasADi_test_Speedy_P.SFunction_P2';
pt(3).baseaddr   = '&CasADi_test_Speedy_P.SFunction_P2';
pt(3).dtname     = 'real_T';


function len = getlenPT
len = 3;

