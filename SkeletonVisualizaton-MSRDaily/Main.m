function Main()
%-----------------------------------------------------------------
%generate mesh
%Author: Jin Qi
%Date:   1/30/2014
%Email:  jqichina@hotmail.com
%copyright2014@gru
%-----------------------------------------------------------------
% subset definition
AS1=[2 3 5 6 10 13 18 20];AS2=[1 4 7 8 9 11 14 12];AS3=[6 14 15 16 17 18 19 20];
AS=AS3;
ASList={AS1 AS2 AS3};
bICA=false;% ICA algorithm
bSparseFeature=true;
%-----------------------------------------------------------------
 %for i=1:numel(ASList)
    % AS=ASList{i}
    Main_Skeleton_MSRAction3D()
    if bICA
       Main_ICA_MSRAction3D();
       Main_HistFeature_ICA_MSRAction3D();
    else
       Main_SparseCoding_Dictionary();
       if bSparseFeature
           Main_HistFeature_SC_MSRAction3D();
       else
           Main_HistFeature_DirectProjection_MSRAction3D();
       end           
    end
    conf=Main_TranTest_MSRAction3D();
%end