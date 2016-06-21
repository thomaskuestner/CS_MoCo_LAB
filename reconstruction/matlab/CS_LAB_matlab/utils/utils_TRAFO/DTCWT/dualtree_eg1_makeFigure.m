
clear
close all


% Set defaults
font_size = 11;
set(0, 'DefaultAxesFontSize',font_size)

line_width = 1;
set(0, 'DefaultLineLineWidth',line_width)

font_size = 15;
set(0, 'DefaultTextFontSize',font_size)




subplot(2,1,1)

dualtree_eg1

box off

% remove defaults
set(0, 'DefaultAxesFontSize','remove')
set(0, 'DefaultLineLineWidth','remove')
set(0, 'DefaultTextFontSize','remove')

if 1
    % paperposition
    paperpos = [2 2 7 5];    % [left, bottom, width, height]

    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', paperpos);
end

print -depsc figures/dualtree_eg1_figure

% then open eps file and save as png file.
