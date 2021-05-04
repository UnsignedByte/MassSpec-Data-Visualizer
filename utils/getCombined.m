function FinalFileOut = getCombined(datasets, datasetnames, parsedProteins, UniqueColumns, UniqueCombineFunctions, UniqueClassFunctions, SingleColumns, SingleClassFunctions) 
    NumFilesRead = length(datasets);

    %Building the main data structure here. All data is read in and retained.
    %This should make it flexible to do alternative readouts on the fly. Should
    %be scalable to any number of input files.
    for i = 1 : NumFilesRead
        TempStruct(i).NameForFile = datasetnames{i};
        TempStruct(i).dat = datasets{i};
    end

    ProteinNamesMap = containers.Map(parsedProteins.fullname, 1:size(parsedProteins,1)); % Create a map from name -> ID (arbitrary #)
    ProteinRanks = [1:size(parsedProteins,1)]; % List of ranks
    ProteinExists = zeros(1,size(parsedProteins,1));

    % disp(TempStruct(i).dat)

    %Match ids if proteins are in groups
    for i = 1:NumFilesRead % loop each file
        jstart = 1; %Save start
        for jreal = 2:size(TempStruct(i).dat,1)+1 % loop through all ranks/proteins in dataset
            if jreal > size(TempStruct(i).dat,1) || TempStruct(i).dat.ProteinRank(jstart) ~= TempStruct(i).dat.ProteinRank(jreal) %if rank is different than before or we are at the end
                %find ids of a "class" of proteins (proteins with the same rank)
                ids = cellfun(@(x) ProteinRanks(ProteinNamesMap(x)), TempStruct(i).dat.Description(jstart:jreal-1));
                ProteinExists(ids) = 1; %set these to exists
                ids = arrayfun(@(x) find(ProteinRanks==x), ids, 'UniformOutput', false);
                ids = [ids{:}];
                ProteinRanks(ids) = min(ProteinRanks(ids)); %take the smallest id to save it to
                jstart = jreal;
            end
        end
    end
    
    % assign new ranks to names list
    % disp(length(find(ProteinExists)))
    % disp(ProteinNames(find(ProteinExists==0)));
    parsedProteins.ranks = ProteinRanks';
    % ProteinNames(:,2) = num2cell(ProteinRanks);
    parsedProteins = parsedProteins(find(ProteinExists), :);

    % sort by new ranks
    parsedProteins = sortrows(parsedProteins,"ranks"); %sort names by id
    parsedProteinsMap = containers.Map(parsedProteins.fullname, 1:size(parsedProteins,1));

    % disp(ProteinNames)

    totUniqueFuncs = zeros(length(UniqueColumns),1);
    totUniqueFuncs(1) = length(UniqueCombineFunctions{1});
    for i = 2:length(UniqueColumns)
        totUniqueFuncs(i) = totUniqueFuncs(i-1)+length(UniqueCombineFunctions{i});
    end

    Result = NaN(size(parsedProteins,1),NumFilesRead*length(UniqueColumns)+totUniqueFuncs(end)+length(SingleColumns)+2);
    ClassResult = [];

    function res = getContaminated(data)
        res = max(2*any(data==2), all(data==1)); %returns 2 if any contaminated (2) exist, or 1 if all are reverse
    end


    realI = 1;
    j = 1;
    for i = 1:parsedProteins.ranks(end)
        startJ = j;
        while j <= size(parsedProteins, 1) && parsedProteins.ranks(j) == i
            for kk = 1:NumFilesRead
                ind = find(strcmp(TempStruct(kk).dat.Description, parsedProteins.fullname(j))); %find index of peptide
                if ~isempty(ind)
                    for k = 1:length(UniqueColumns)
                        Result(j,(k-1)*NumFilesRead+totUniqueFuncs(k)-length(UniqueCombineFunctions{k})+kk) = max(table2array(TempStruct(kk).dat(ind,UniqueColumns(k))), [], 1); %paste data into result
                    end
                    %unique columns should be equal across all datasets, so only
                    %run for all
                    for k = 1:length(SingleColumns)
                        Result(j,length(UniqueColumns)*NumFilesRead+totUniqueFuncs(end)+k) = max(max(Result(j,length(UniqueColumns)*NumFilesRead+totUniqueFuncs(end)+k), TempStruct(kk).dat{ind, SingleColumns(k)}), [], 1); %paste data into result
                    end
                end
            end
            % add combined columns
            for k = 1:length(UniqueColumns)
                currSliceInd = k*NumFilesRead+totUniqueFuncs(k)-length(UniqueCombineFunctions{k});
                for l = 1:length(UniqueCombineFunctions{k})
                    Result(j, currSliceInd+l)=feval(UniqueCombineFunctions{k}{l}, Result(j, currSliceInd-NumFilesRead+1:currSliceInd));
                end
            end

            % code to find the reverses and contaminants and flag them
            %for removal.
            Result(j, end-1) = 0;
            % If it has any contams, skip
            if       ~(isempty(strfind(parsedProteins.fullname(j),'Common contaminant')) ...
                    && isempty(strfind(parsedProteins.fullname(j),'eratin, type')) ...
                    && isempty(strfind(parsedProteins.fullname(j),' desmo')) ...
                    && isempty(strfind(parsedProteins.fullname(j),'dermi')) ...
                    && isempty(strfind(parsedProteins.fullname(j),'plak')))
                Result(j, end-1) = 2;
            elseif ~ismissing(parsedProteins.reverse(j))
                % Contam iff no >tr and no >sp
                Result(j, end-1) = 1;
            end
            Result(j, end) = i;
            j = j + 1;
        end
        j = j - 1;
        if j >= startJ
            ClassResult(end+1,:) = NaN(1, size(Result,2));
            for k = 1:length(UniqueColumns)
                currSliceInd = (k-1)*NumFilesRead+totUniqueFuncs(k)-length(UniqueCombineFunctions{k});
                for kk = 1:NumFilesRead
                    ClassResult(end, currSliceInd+kk) = feval(UniqueClassFunctions{k}, Result(startJ:j,currSliceInd+kk));
                end
                Slice = ClassResult(end,currSliceInd+1:currSliceInd+NumFilesRead);
                for l = 1:length(UniqueCombineFunctions{k})
                    ClassResult(end, k*NumFilesRead+totUniqueFuncs(k)-length(UniqueCombineFunctions{k})+l)=feval(UniqueCombineFunctions{k}{l}, Slice);
                end
            end
            for k = 1:length(SingleColumns)
                currInd = length(UniqueColumns)*NumFilesRead+totUniqueFuncs(end)+k;
                ClassResult(end,currInd) = feval(SingleClassFunctions{k}, Result(startJ:j,currInd));
            end
            ClassResult(end,end-1) = getContaminated(Result(startJ:j, end-1)); %add flags
            ClassResult(end,end) = i;
            realI = realI + 1;
        end
        j = j + 1;
    end

    sortOrd = [size(ClassResult,2), -NumFilesRead-4, -NumFilesRead-5];
    ClassResult = sortrows(ClassResult, [size(ClassResult,2)-1, -NumFilesRead-1, -NumFilesRead-2]);

    Result = num2cell(Result);
    Result = [Result(:,end) cellstr(parsedProteins.fullname) Result(:,1:end-1)];
    ClassResult = num2cell(ClassResult);
    ClassResult = [ClassResult(:,end) cell(size(ClassResult,1),2) ClassResult(:,1:end-1)];

    CombinedRes = cell(size(ClassResult,1)+size(Result,1), size(ClassResult,2)+1);
    realI = 1;
    
    for i = 1:size(ClassResult,1)
        CombinedRes(realI,:) = [i ClassResult(i,2:end) 1];
        found = Result(cell2mat(Result(:,1))==ClassResult{i,1},:); %get rows in class
        fnames = found(:,2); %take out names
        pranks =arrayfun(@(x) parsedProteinsMap(x), string(fnames));

        found = [found(:,2) cellstr(parsedProteins.genename(pranks)) found(:,3:end)];
        found(:,1) = num2cell(1:size(found,1));
        found = sortrows(found, sortOrd); %sort by rows
        CombinedRes{realI,2} = char(strjoin(unique(parsedProteins.proteinname(pranks(~cellfun('isempty',parsedProteins.proteinname(pranks))))), '/')); %name of class
        CombinedRes{realI,3} = char(strjoin(unique(parsedProteins.proteinname(pranks(~cellfun('isempty',parsedProteins.genename(pranks))))), '/')); %name of class
        % disp(found)
        found(:,1) = fnames(cell2mat(found(:,1))); %replace names
        CombinedRes(realI+1:realI+size(found,1),2:end-1) = found;
        CombinedRes(realI+1:realI+size(found,1),1) = {i};
        CombinedRes(realI+1:realI+size(found,1),end) = {0};
        realI = realI+size(found,1)+1;
    end
    %currently sorts by max then sum of spectra

    % If this search is for the Jackson Lab, use the following segment:

    %     for j = 1:length(craptxt(:,1))
    %         if isempty(strcmp(FileOutCombined{i,1},craptxt{j,1})) == 0
    %             reversefinder(i,1) = 1;
    %         end;
    %     end;
    % end

    % FileOutFinalClasses = cell(FileOutFinal{end,1},size(FileOutFinal,2)); %get header classes
    % 
    % j = 1;
    % for i = 1:size(FileOutFinalClasses,1)
    %     chunk = {};
    %     while FileOutFinal{j,1} == i
    %         chunk(end+1, :) = FileOutFinal(j,:);
    %         j = j + 1;
    %     end
    %     class = cell(1,size(chunk,2));
    %     for ll = 1:length(UniqueClassFunctions)
    %         for k = 1:NumFilesRead
    %             
    %         end
    %     end
    %     for k = 2:size(chunk,2)
    %         class{k} = UniqueClassFunctions{
    %     end
    % end

    FileOutFinalDeComma = CombinedRes;

    for i = 1:size(CombinedRes, 1)
        if isempty(strfind(CombinedRes{i,2},',')) == 0
            FileOutFinalDeComma{i,2} = strrep(CombinedRes{i,2}, ',', '-');
        end
    end

    %Here, I add in the file names in a readable format so they can be column
    %headers. Note that these column headers are a tad different than before so
    %they are readable, and should be based on the input file names as much as
    %possible.

    HeaderFileString = cell(size(CombinedRes,2), 1);
    HeaderFileString(1:3) = {'Rank_Number'; 'Protein_Name'; 'Gene_Name'};
    HeaderFileString{end-1} = 'Contaminant';
    HeaderFileString{end} = 'Row_Type';

    for i = 1 : NumFilesRead
%         CleanedFileName = strrep(CleanedFileName, '_control',''); %Remove after I'm done with Janos' dataset
%         CleanedFileName = strrep(CleanedFileName, '25fmol_ul_6x5spike_',''); %Remove after I'm done with Janos' dataset
        %label by dataset
        for j = 1:length(UniqueColumns)
            tval = [TempStruct(1).dat.Properties.VariableNames{UniqueColumns(j)} '_' TempStruct(i).NameForFile];
            HeaderFileString{3+(j-1)*NumFilesRead+totUniqueFuncs(j)-length(UniqueCombineFunctions{j})+i} = matlab.lang.makeValidName(tval);
        end
    end

    %label single columns
    for j = 1:length(SingleColumns)
        HeaderFileString{3+length(UniqueColumns)*NumFilesRead+totUniqueFuncs(end)+j} = TempStruct(1).dat.Properties.VariableNames{SingleColumns(j)};
    end

    %label max and sum
    for j = 1:length(UniqueColumns)
        for jj = 1:length(UniqueCombineFunctions{j})
            HeaderFileString{3+j*NumFilesRead+totUniqueFuncs(j)-length(UniqueCombineFunctions{j})+jj} = [UniqueCombineFunctions{j}{jj} '_' TempStruct(1).dat.Properties.VariableNames{UniqueColumns(j)}];
        end
    end

    %Now we make the final output file for printing to a csv. we first make a
    %table rearrange it appropriately, and then have it fill in the rank number
    %values before printing_ it out.

    FinalFileOut = cell2table(FileOutFinalDeComma);
    FinalFileOut.Properties.VariableNames = HeaderFileString;

    
%     fid = fopen([filename '.json'], 'w');
%     fprintf(fid, jsonencode(FinalFileOut,'ConvertInfAndNaN', false));
%     fclose(fid);
%     output = table2struct(FinalFileOut);
%     save('output.mat', 'output', '-v7.3');
end