function lc_save_png(figH, fname)
    if exist('exportgraphics', 'file') == 2
        exportgraphics(figH, fname, 'Resolution', 300);
    else
        print(figH, fname, '-dpng', '-r200');
    end
end
