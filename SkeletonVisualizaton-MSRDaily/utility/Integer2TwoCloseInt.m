
function [n1 n2]=Integer2TwoCloseInt(n)
%------------------------------------------------------
%Author: Jin Qi
%Date:   11/24/2014
%Email:  jqichina@hotmail.com
%copyright2014@cnmc
%------------------------------------------------------
FactorNum=factor(n);
% n1=1;n2=1;
DeltaDiff=n; % difference between two factors
if length(FactorNum)==1
    n1=1;n2=FactorNum;
else
    for i=1:length(FactorNum)-1
        n1=prod(FactorNum(1:i));
        n2=prod(FactorNum(i+1:end));
        if abs(n1-n2)<DeltaDiff
            DeltaDiff=abs(n1-n2);
        else
            n1=prod(FactorNum(1:i-1));
            n2=prod(FactorNum(i:end));
        end
    end
end
end