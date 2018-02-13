function update_network(A,S_var,h)

[L M]=size(A);


figure(1)

dim=ceil(sqrt(M));
[dim1 dim2]=Integer2TwoCloseInt(M)
% sz=sqrt(L);
[sz1 sz2]=Integer2TwoCloseInt(L)

for i=1:M
  a_max=max(abs(A(:,i)));
  set(subplot(dim1,dim2,i),'CLim',[-a_max a_max])
  set(h(i),'CData',reshape(A(:,i),sz1,sz2))
  drawnow
end


if exist('S_var','var')
  figure(2)
  subplot(211), bar(S_var), title('s variance')
  subplot(212), bar(sqrt(sum(A.*A))), title('basis norm (L2)')
end

drawnow
