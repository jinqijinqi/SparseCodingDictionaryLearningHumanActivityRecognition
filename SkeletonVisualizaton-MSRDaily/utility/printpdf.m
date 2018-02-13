function printpdf(h,outfilename)
set(h,'units','normalized','outerposition',[0 0 1 1]);
set(h, 'PaperSize', [6.25 7.5]);
set(h, 'PaperPositionMode', 'manual');
set(h, 'PaperPosition', [0 0 6.25 7.5]);

set(h, 'PaperUnits', 'inches');
set(h, 'PaperSize', [6.25 7.5]);
set(h, 'PaperPositionMode', 'manual');
set(h, 'PaperPosition', [0 0 6.25 7.5]);


set(h, 'renderer', 'painters');
print(h, '-dpdf', 'outfilename.pdf');
print(h, '-dpng', 'my-figure.png');
print(h, '-depsc2', 'my-figure.eps');

print(h, '-dpdf', [outfilename '.pdf']);
print(h, '-dpng', [outfilename '.png']);
print(h, '-depsc2',[outfilename '.eps']);
mycmd=['pdfcrop ' outfilename '.pdf'];
system(mycmd);
% !pdfcrop [outfilename '.pdf'] [outfilename '-crop.pdf'];