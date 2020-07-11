function params = getLatestParams(p)
  % datetime format
	fmt = 'dd-MM-yyyy_HH:mm:ss';

	list = dir([p, '/*.txt']);
	% cell array of the times
	parsedTokens = regexp({list.name}, '(.+)_(\d\d-\d\d-\d{4}_\d\d:\d\d:\d\d)\.txt', 'tokens'); % list of parsed tokens
	parsedTokens = horzcat(parsedTokens{:}); % unnest inner arrays
	parsedTokens = vertcat(parsedTokens{:}); % list of 1x2 cells to 10x2 matrix
	parsedTokens(:,2) = num2cell(cellfun(@(x) datetime(x, 'InputFormat', fmt), parsedTokens(:,2))); % convert objects to datetime

	% get max datetime (latest) for each type of param file
	params = cellfun(@(x) {x, [x '_' datestr(max([parsedTokens{contains(parsedTokens(:,1), x),2}]), fmt) '.txt']}, unique(parsedTokens(:,1)), 'UniformOutput', false);
	params = vertcat(params{:}); % unnest cells again
end