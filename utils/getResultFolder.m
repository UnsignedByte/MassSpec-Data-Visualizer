function folderName = getResultFolder(filename)
    regStr = '^\d+_(.+?)_(.+?)_';
    toks = regexp(filename, regStr, 'tokens');
    folderName = [toks{:}{1} '_' toks{:}{2}];
    if ~isfolder('Results')
        mkdir('Results')
    end
    if ~isfolder(fullfile('Results', folderName))
        mkdir(fullfile('Results', folderName))
    end
end