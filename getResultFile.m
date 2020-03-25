function filename = getResultFile(filename)
    regStr = '^\d+_.+?_.+?_.+?_(.+)\.raw';
    toks = regexp(filename, regStr, 'tokens');
    filename = toks{:}{1};
    
end