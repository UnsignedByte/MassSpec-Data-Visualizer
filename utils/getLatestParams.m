function params = getLatestParams(p)
  % datetime format

	list = dir([p, '/*.txt']);
	% cell array of the times
	parsedTokens = regexp({list.name}, '(.+)_((?:\d\d-){2}\d{4}_(?:\d\d:){2}\d\d)\.txt', 'tokens'); % list of parsed tokens
	parsedTokens = horzcat(parsedTokens{:}); % unnest inner arrays
	parsedTokens = vertcat(parsedTokens{:}); % list of 1x2 cells to 10x2 matrix
	parsedTokens(:,2) = num2cell(cellfun(@(x) datetime(x, 'InputFormat', 'dd-MM-yyyy_HH:mm:ss'), parsedTokens(:,2))); % convert objects to datetime

	% return map containing filename for the latest version of each param
	keyset = unique(parsedTokens(:,1));
	params = containers.Map(keyset, cellfun(@(x) [x '_' datestr(max([parsedTokens{contains(parsedTokens(:,1), x),2}]), 'dd-mm-yyyy_HH:MM:SS') '.txt'], keyset, 'UniformOutput', false));
end