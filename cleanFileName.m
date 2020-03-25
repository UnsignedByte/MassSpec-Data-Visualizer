function CleanedFileName = cleanFileName(fname)
    CleanedFileName = strrep(fname, '.raw_20', '_');
    CleanedFileName = strrep(CleanedFileName, '_Byonic','');
    CleanedFileName = strrep(CleanedFileName, '_Elite','');
end