function mypdf(fname,r,s)

%disp('MYPDF CURRENTLY DISABLED');
%return


if nargin < 2,
    r = .71; % height/width ratio
end;
if nargin < 3, 
    s = 1; % scaling of font size
end;


set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperSize',s*[13,r*13]);
set(gcf,'PaperPosition',s*[0,0,13,r*13]);
print(gcf,'-dpdf',  fname);
print(gcf,'-depsc2', fname);
print(gcf,'-dpng','-r200', fname);
saveas(gcf,[  fname '.fig'])
%matlab2tikz([fname '.tikz'], 'height', '\figureheight', 'width', '\figurewidth' )