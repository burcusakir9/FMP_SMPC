function saveFigure(f, cfg, name)
%SAVEFIGURE  Save figure f as figures/<name>.png (200 dpi).
% Retries once if the file is locked (e.g. open in a viewer), then falls
% back to <name>_new.png so a run never fails on a locked image.
    file = fullfile(cfg.dirFigures, [name '.png']);
    for attempt = 1:2
        try
            exportgraphics(f, file, 'Resolution', 200);
            return;
        catch
            pause(2);
        end
    end
    alt = fullfile(cfg.dirFigures, [name '_new.png']);
    warning('saveFigure: %s is locked, saved to %s instead.', file, alt);
    exportgraphics(f, alt, 'Resolution', 200);
end
