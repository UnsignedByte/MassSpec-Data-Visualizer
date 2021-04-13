function parsed = parseProteins(proteins, fmts) %return info about each protein in a list of protein names
    if ~exist('fmts', 'var')
        fmts = {};
    end
    fmts = [fmts,...
            {
            '^{>(?<reverse>Reverse)\s*}{>(?<tr>tr)\s*}{>(?<sp>sp)}\|(?<dbname>.+?)\|(?<proteinname>.+?){_(?<organism>.+?)}\s(?<type>.+?){\sOS=(?<os>(?:.|\s|)+?)}{\sOX=(?<ox>(?:.|\s|)+?)}{\sGN=(?<genename>(?:.|\s|)+?)}{\sPE=(?<pe>(?:.|\s|)+?)}{\sSV=(?<sv>(?:.|\s|)+?)}$',...
            '^>{()}{\s(?<proteinname>(.|\s)+)\s}{\s\[(?<organism>(?:.|\s)+)\]}$',...
            '^>?(?<reverse>Reverse)\s*>?(?<tr>tr)\s*>?(?<sp>sp)(?:\||\s)(?<proteinname>.+?)$'
            }];

    function m = gparse(x)
        rng(0, 'twister'); % reset rng seed so we can use this later
        m = regexprep(x, '{(.+?)}', '(?<generated_${num2str(floor(rand()*1e6))}>$1)?');
    end

    fmtsR = cellfun(@(x) gparse(x), fmts, 'UniformOutput', false); % Convert custom regex to usable regex
    fmtsGroups = cellfun(@(x) regexp(x, '{(.+?)}', 'tokens'), fmts, 'UniformOutput', false); % Save groups inside {} for later
    
    %Identifier, Gene name, 
    %disp(proteins)
    function p = parse(x)
        p = struct;
        for i = 1:length(fmts)
            if regexp(x, fmtsR{i})
                p = regexp(x, fmtsR{i}, 'names');
                rng(0, 'twister');
                for j = 1:length(fmtsGroups{i})
                    fld = ['generated_' num2str(floor(rand()*1e6))];
                    p = mergeStruct(p, regexp(p.(fld), fmtsGroups{i}{j}{1}, 'names'));
                    p = rmfield(p, fld);
                end
                p.fullname = x;
                return;
            end
        end
        disp("NO MATCH FOUND")
    end

    parsed = cellfun(@(x) parse(x), proteins, 'UniformOutput', false);
    % parsed(:, 9:end) = cellfun(@(x) x(5:end), parsed(:, 9:end), 'UniformOutput', false); %remove SV=, etc 
end