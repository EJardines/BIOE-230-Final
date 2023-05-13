
function [] = PlotSettingsLogLog(fontsize,fh)

    set(fh,'color','white');
    box off % turns right and top portions of the box off (leaves axes)
    A = gca;
%     set(fh,'color','none');
    A.FontSize = fontsize;
    A.TickLabelInterpreter = 'latex';
    A.XScale = 'log';
    A.YScale = 'log';
%     set(gca, 'XScale', 'log', 'YScale', 'log');
%     fprintf('Successfully ran the plotting function!');
    
end