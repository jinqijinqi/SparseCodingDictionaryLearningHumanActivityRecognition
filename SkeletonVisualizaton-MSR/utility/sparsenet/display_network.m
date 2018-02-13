function h=display_network(A,S_var)

[L M]=size(A);


figure(1)

% dim=ceil(sqrt(M));
[dim1 dim2]=Integer2TwoCloseInt(M);
% sz=sqrt(L);
[sz1 sz2]=Integer2TwoCloseInt(L)

h=zeros(M,1);

for i=1:M
  a_max=max(abs(A(:,i)));
  subplot(dim1,dim2,i)
  h(i)=imagesc(reshape(A(:,i),sz1,sz2),[-a_max a_max]);
%   h(i)=imagesc(reshape(A(:,i),sz1,sz2),'EraseMode','none',[-a_max a_max]);
  axis square, axis off
  drawnow
end


if exist('S_var','var')
  figure(2)
  subplot(211), bar(S_var), title('s variance')
  subplot(212), bar(sqrt(sum(A.*A))), title('basis norm (L2)')
end

drawnow
end
% function [n1 n2]=Integer2TwoCloseInt(n)
% FactorNum=factor(n);
% % n1=1;n2=1;
% DeltaDiff=n; % difference between two factors
% if length(FactorNum)==1
%     n1=1;n2=FactorNum;
% else
%     for i=1:length(FactorNum)-1
%         n1=prod(FactorNum(1:i));
%         n2=prod(FactorNum(i+1:end));
%         if abs(n1-n2)<DeltaDiff
%             DeltaDiff=abs(n1-n2);
%         else
%             n1=prod(FactorNum(1:i-1));
%             n2=prod(FactorNum(i:end));
%         end
%     end
% end
% end
