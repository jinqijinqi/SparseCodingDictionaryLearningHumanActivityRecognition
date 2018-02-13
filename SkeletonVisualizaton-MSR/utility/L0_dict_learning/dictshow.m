function F = dictshow(D)
%% Show columns of dictionary matrix D as blocks
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Use of this code is free for research purposes only.
%
%Author:  Yuhui Quan
%
%Last Revision: 22-Jun-2014
 %
 %
    p = 3.5;
    
    [dim_atom, n_atom] = size(D);
    sz_edge = sqrt(dim_atom);
    min_D = min(D(:));
    if (min_D >= 0)
       m = 0;
       sig = sqrt(mean(D(:).^2));
    else
       m = mean(D(:));
       sig = sqrt(mean((D(:)-m).^2));
    end
    D = D - m;
    D = min(max(D,-p*sig),p*sig);
    max_D = max((D(:)));
    min_D = min((D(:)));
    D = (D - min_D)/(max_D - min_D);
   
    n_rbin = floor(sqrt(n_atom));
    n_cbin = n_rbin;
    
    F = zeros((sz_edge+1)*n_rbin+1,(sz_edge+1)*n_cbin+1); %Fig
    for i = 1: n_rbin
        for j = 1: n_cbin
          P = reshape(D(1:dim_atom,(i-1)*n_cbin+j), [sz_edge,sz_edge]);
          F((i-1)*(sz_edge+1)+2:i*(sz_edge+1),...
             (j-1)*(sz_edge+1)+2:j*(sz_edge+1),:)=P;
        end
    end

    colormap('bone');
    imagesc(F);

end

