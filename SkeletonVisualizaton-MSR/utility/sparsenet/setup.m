%% Setup Matlab according to Bruno's notes
%% CHA Feb. 27, 2002
load IMAGES;
A = rand(64)-0.5;
A = A*diag(1./sqrt(sum(A.*A)));
figure(1); colormap(gray);
figure(2);