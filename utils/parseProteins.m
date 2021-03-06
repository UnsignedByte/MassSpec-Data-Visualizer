function parsed = parseProteins(proteins) %return info about each protein in a list of protein names
    regStr = '^>?(.+?)\|(.+?)\|(.+?)_(.+?)\s((?:.+?)+?)(\sOS=(?:.|\s|)+?)(\sOX=\d+?)(\sGN=.+?)?(\sPE=\d+)?(\sSV=\d+)?$';

    %disp(proteins)
    
    function p = parse(x)
    	if regexp(x, regStr)
    		p = regexp(x, regStr, 'tokens');
    		p{1}(1,6:10) = cellfun(@(x) split([x '='], '='), p{1}(1, 6:10), 'UniformOutput', false); %Split by = sign
    		p{1}(1,6:10) = cellfun(@(x) x{2}, p{1}(1, 6:10), 'UniformOutput', false); %Take value only
    	else
    		p = cell(1,10);
    		p{3} = regexp(x, '^>?(.+?)$', 'tokens');
    		p{3} = p{3}{1}{1};
    		p = {p};
    	end
    end

    parsed = cellfun(@(x) [x, parse(x), {cell(1,10)}], proteins, 'UniformOutput', false);
    parsed = cellfun(@(x) [x{1}, x{2}], parsed, 'UniformOutput', false);
    parsed = vertcat(parsed{:});
    % parsed(:, 9:end) = cellfun(@(x) x(5:end), parsed(:, 9:end), 'UniformOutput', false); %remove SV=, etc 
end