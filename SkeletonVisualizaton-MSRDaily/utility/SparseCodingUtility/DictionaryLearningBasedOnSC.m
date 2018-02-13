function [LearnedDict WhitenMatrix]=DictionaryLearningBasedOnSC(NumAtom,TrainData)
%------------------------------------------------------
%Author: Jin Qi
%Date:   11/24/2014
%Email:  jqichina@hotmail.com
%copyright2014@cnmc
%------------------------------------------------------
%% 
% learn a dictioary based on sparse coding
% input:
%         NumAtom? number of atoms in dictionary
%         TrainData: training data
% output:
%         Dict: learned dictionary
%% initialize dictionary
X=TrainData;
[WhitenMatrix WhiteData]=PCAWhitenData(X);% whiten data
X=WhiteData;
WhiteLength=size(WhiteData,1); % length of sample becomes smaller after whiten
A = rand(WhiteLength,NumAtom)-0.5;
A = A*diag(1./sqrt(sum(A.*A)));
%% 



num_trials=10;
batch_size=size(X,2);
% batch_size=100;

% num_images=size(IMAGES,2);
% image_size=sqrt(size(IMAGES,1));
% BUFF=4;

[L M]=size(A);
sz=sqrt(L);

eta = 1.0; % original
% eta = 0.1;
noise_var= 0.01;%orginal
% noise_var=0.001;
beta= 2.2;%orginal
% beta=0.2
sigma=0.316; %orginal
% sigma=0.03
tol=.01; %orginal
% tol=0.001;

VAR_GOAL=0.1;
S_var=VAR_GOAL*ones(M,1);
var_eta=.001;
alpha=.02;
gain=sqrt(sum(A.*A))';

% X=zeros(L,batch_size);

display_every=5;
NAtom=20;% display number of atoms
if (exist('disp_handle','var'))
  update_network(A,S_var,disp_handle);
else
  disp_handle=display_network(A(:,1:NAtom),S_var);
end

for t=1:num_trials
    fprintf('Dictionary learning: %d(%d)\n',t,num_trials);

% choose an image for this batch

%   i=ceil(num_images*rand);
%   this_image=reshape(IMAGES(:,i),image_size,image_size)';
%   
% % extract subimages at random from this image to make data vector X
% 
%   for i=1:batch_size
%     r=BUFF+ceil((image_size-sz-2*BUFF)*rand);
%     c=BUFF+ceil((image_size-sz-2*BUFF)*rand);
%     X(:,i)=reshape(this_image(r:r+sz-1,c:c+sz-1),L,1);
%   end

% calculate coefficients for these data via conjugate gradient routine

  S=cgf_fitS(A,X,noise_var,beta,sigma,tol);

% calculate residual error

  E=X-A*S;

% update bases

  dA=zeros(L,M);
  for i=1:batch_size
    dA = dA + E(:,i)*S(:,i)';
  end
  dA = dA/batch_size;

  A = A + eta*dA;

% normalize bases to match desired output variance

  for i=1:batch_size
    S_var = (1-var_eta)*S_var + var_eta*S(:,i).*S(:,i);
  end
  gain = gain .* ((S_var/VAR_GOAL).^alpha);
  normA=sqrt(sum(A.*A));
  for i=1:M
    A(:,i)=gain(i)*A(:,i)/normA(i);
  end

% display

  if (mod(t,display_every)==0)
%     update_network(A(:,1:NAtom),S_var,disp_handle);
%     disp_handle=display_network(A(:,1:NAtom),S_var);
  end
end
LearnedDict=A';
