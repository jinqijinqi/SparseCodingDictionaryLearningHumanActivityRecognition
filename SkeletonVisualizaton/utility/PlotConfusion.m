function PlotConfusion(Confusion,TestLabel,PredictLabel,FontSize)
%generate mesh
%Author: Jin Qi
%Date:   4/14/2014
%Email:  jqichina@hotmail.com
%copyright2014@gru
%%
FS=FontSize;
mat=Confusion;

% mat = rand(5);           %# A 5-by-5 matrix of random values from 0 to 1
% mat(3,3) = 0;            %# To illustrate
% mat(5,2) = 0;            %# To illustrate
imagesc(mat);            %# Create a colored plot of the matrix values
colormap(flipud(gray));  %# Change the colormap to gray (so higher values are
                         %#   black and lower values are white)

textStrings = num2str(mat(:),'%1.f');  %# Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding

%% ## New code: ###
idx = find(strcmp(textStrings(:), '0'));
textStrings(idx) = {'   '};
%% ################
[Row Col]=size(mat);
[x,y]=meshgrid(1:Col,1:Row);
% [x,y] = meshgrid(1:5);   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                'HorizontalAlignment','center','fontsize',FS);
% set(hStrings,'FontSize',16);
midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
textColors = repmat(mat(:) > midValue,1,3);  %# Choose white or black for the
                                             %#   text color of the strings so
                                             %#   they can be easily seen over
                                             %#   the background color
set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
set(gca,'XTick',1:Col,...                         %# Change the axes tick marks
        'XTickLabel',TestLabel,...  %#   and tick labels
        'YTick',1:Row,...
        'YTickLabel',PredictLabel,...
        'TickLength',[0 0],'fontsize',FS);
v=get(gca,'Position'); % set axis
% set(gca,'Position',[v(1)*1 v(2:4)]);
set(gca,'Position',[v(1)*1.1 v(2) v(3)*1.5 v(4)]);
xticklabel_rotate([1:numel(PredictLabel)],45,PredictLabel,'interpreter','none');
% set(gca,'XTick',1:5,...                         %# Change the axes tick marks
%         'XTickLabel',{'A','B','C','D','E'},...  %#   and tick labels
%         'YTick',1:5,...
%         'YTickLabel',{'A','B','C','D','E'},...
%         'TickLength',[0 0]);