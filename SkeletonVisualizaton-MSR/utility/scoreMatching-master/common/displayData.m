function [h, display_array] = displayData(X, example_width)
%DISPLAYDATA Display 2D data in a nice grid
%   [h, display_array] = DISPLAYDATA(X, example_width) displays 2D data
%   stored in X in a nice grid. It returns the figure handle h and the 
%   displayed array if requested. Expects examples in X to be in rows.

% Set example_width automatically if not passed in
if ~exist('example_width', 'var') || isempty(example_width) 
	example_width = round(sqrt(size(X, 2)));
end

% Gray Image
colormap(gray);

% Compute rows, cols
[m n] = size(X);
% fs=factor(size(X,2));
% fd=factor(size(X,1)); % factor number of dictionary
[example_width example_height]=Factor2Int(n);
[display_rows display_cols]=Factor2Int(m);
% if length(fs)==1
%     example_width=1;
%     example_height =(n / example_width);
%     display_rows=1;
%     display_cols=m/display_rows;
% else
%     example_width=prod(fs(1:end-1));
%     example_height =(n / example_width);
%     display_rows=prod(fd(1:end-1));
%     display_cols=m/display_rows;
% end

% Compute number of items to display
% display_rows = floor(sqrt(m));
% display_cols = ceil(m / display_rows);

% Between images padding
pad = 1;

% Setup blank display
display_array = - ones(pad + display_rows * (example_height + pad), ...
                       pad + display_cols * (example_width + pad));

% Copy each example into a patch on the display array
curr_ex = 1;
for j = 1:display_rows
	for i = 1:display_cols
% 		if curr_ex > m, 
        if curr_ex > n,
			break; 
		end
		% Copy the patch
		
		% Get the max value of the patch
		max_val = max(abs(X(curr_ex, :)));
		display_array(pad + (j - 1) * (example_height + pad) + (1:example_height), ...
		              pad + (i - 1) * (example_width + pad) + (1:example_width)) = ...
						reshape(X(curr_ex, :), example_height, example_width) / max_val;
		curr_ex = curr_ex + 1;
	end
% 	if curr_ex > m, 
    if curr_ex > n,
		break; 
	end
end

% Display Image
h = imagesc(display_array, [-1 1]);

% Do not show axis
axis image off

drawnow;

end
function [a b]=Factor2Int(x)
sr=floor(sqrt(x));
while true
    if mod(x,sr)==0
       a=sr;b=x/sr;
       return;
    else
        sr=sr-1;
    end
end
end
