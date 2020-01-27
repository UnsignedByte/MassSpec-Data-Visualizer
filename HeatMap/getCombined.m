function getCombined(datasets, datasetnames, UniqueColumns, UniqueCombineFunctions, UniqueClassFunctions, SingleColumns, SingleClassFunctions, filename) 
    NumFilesRead = length(datasets);
    ProteinNames = []; %all protein names

    %Building the main data structure here. All data is read in and retained.
    %This should make it flexible to do alternative readouts on the fly. Should
    %be scalable to any number of input files.
    for i = 1 : NumFilesRead
        TempStruct(i).NameForFile = datasetnames{i};
        dat = datasets{i};
        dat(:,3:end) = fillmissing(dat(:,3:end), 'constant', 0); %replace NaN with zero
        TempStruct(i).dat = dat;
        ProteinNames = [ProteinNames; dat.Description];
    end

    ProteinNames = unique(ProteinNames); %get only unique
    ProteinNamesMap = containers.Map(ProteinNames, 1:length(ProteinNames));
    ProteinRanks = [1:length(ProteinNames)];

    %Match ids if proteins are in groups
    for i = 1:NumFilesRead
        j = 1;
        jstart = 1;
        for jreal = 2:size(TempStruct(i).dat,1) %loop through all ranks in dataset
            j = TempStruct(i).dat.ProteinRank(jreal);
            if TempStruct(i).dat.ProteinRank(jstart) ~= j %if rank changed
                %find ids of a "class" of proteins (proteins with the same rank)
                ids = cell2mat(cellfun(@(x) ProteinNamesMap(x), TempStruct(i).dat.Description(jstart:jreal-1), 'UniformOutput', false));
                ProteinRanks(ids) = min(ProteinRanks(ids));
                jstart = jreal;
            end
        end
    end
    
    ProteinNames(:,2) = num2cell(ProteinRanks); 
    
    ProteinNames = sortrows(ProteinNames,2); %sort names by id

    toc;

    totUniqueFuncs = zeros(length(UniqueColumns),1);
    totUniqueFuncs(1) = length(UniqueCombineFunctions{1});
    for i = 2:length(UniqueColumns)
        totUniqueFuncs(i) = totUniqueFuncs(i-1)+length(UniqueCombineFunctions{i});
    end

    Result = NaN(length(ProteinNames),NumFilesRead*length(UniqueColumns)+totUniqueFuncs(end)+length(SingleColumns)+2);
    ClassResult = [];

    realI = 1;
    j = 1;
    for i = 1:ProteinNames{end,2}
        startJ = j;
        while j <= size(ProteinNames, 1) && ProteinNames{j,2} == i
            for kk = 1:NumFilesRead
                ind = find(ismember(TempStruct(kk).dat.Description, ProteinNames{j,1})); %find index of peptide
                if ~isempty(ind)
                    for k = 1:length(UniqueColumns)
                        Result(j,(k-1)*NumFilesRead+totUniqueFuncs(k)-length(UniqueCombineFunctions{k})+kk) = table2array(TempStruct(kk).dat(ind,UniqueColumns(k))); %paste data into result
                    end
                    %unique columns should be equal across all datasets, so only
                    %run for all
                    for k = 1:length(SingleColumns)
                        Result(j,length(UniqueColumns)*NumFilesRead+totUniqueFuncs(end)+k) = max(Result(j,length(UniqueColumns)*NumFilesRead+totUniqueFuncs(end)+k), TempStruct(kk).dat{ind, SingleColumns(k)}); %paste data into result
                    end
                end
            end
            % add combined columns
            for k = 1:length(UniqueColumns)
                currSliceInd = k*NumFilesRead+totUniqueFuncs(k)-length(UniqueCombineFunctions{k});
                for l = 1:length(UniqueCombineFunctions{k})
                    Result(j, currSliceInd+l)=UniqueCombineFunctions{k}{l}(Result(j, currSliceInd-NumFilesRead+1:currSliceInd));
                end
            end

            % code to find the reverses and contaminants and flag them
            %for removal.
            Result(j, end-1) = 0;
            if       ~(isempty(strfind(ProteinNames{j,1},'>Reverse')) ...
                    && isempty(strfind(ProteinNames{j,1},'Common contaminant')) ...
                    && isempty(strfind(ProteinNames{j,1},'eratin, type')) ...
                    && isempty(strfind(ProteinNames{j,1},' desmo')) ...
                    && isempty(strfind(ProteinNames{j,1},'dermi')) ...
                    && isempty(strfind(ProteinNames{j,1},'plak')))
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
                    ClassResult(end, currSliceInd+kk) = UniqueClassFunctions{k}(Result(startJ:j,currSliceInd+kk));
                end
                Slice = ClassResult(end,currSliceInd+1:currSliceInd+NumFilesRead);
                for l = 1:length(UniqueCombineFunctions{k})
                    ClassResult(end, k*NumFilesRead+totUniqueFuncs(k)-length(UniqueCombineFunctions{k})+l)=UniqueCombineFunctions{k}{l}(Slice);
                end
            end
            for k = 1:length(SingleColumns)
                currInd = length(UniqueColumns)*NumFilesRead+totUniqueFuncs(end)+k;
                ClassResult(end,currInd) = SingleClassFunctions{k}(Result(startJ:j,currInd));
            end
            ClassResult(end,end-1) = any(Result(startJ:j, end-1)); %add flags
            ClassResult(end,end) = i;
            realI = realI + 1;
        end
        j = j + 1;
    end

    sortOrd = [size(ClassResult,2), -NumFilesRead-3, -NumFilesRead-4];
    ClassResult = sortrows(ClassResult, [size(ClassResult,2)-1, -NumFilesRead-1, -NumFilesRead-2]);

    Result = num2cell(Result);
    Result = [Result(:,end) ProteinNames(:, 1) Result(:,1:end-1)];
    ClassResult = num2cell(ClassResult);
    ClassResult = [ClassResult(:,end) cell(size(ClassResult,1),1) ClassResult(:,1:end-1)];

    CombinedRes = cell(size(ClassResult,1)+size(Result,1), size(Result,2)+1);
    realI = 1;

    for i = 1:size(ClassResult,1)
        CombinedRes(realI,:) = [i ClassResult(i,2:end) 1];
        found = Result(cell2mat(Result(:,1))==ClassResult{i,1},:); %get rows in class
        fnames = found(:,2); %take out names
        found = found(:,2:end);
        found(:,1) = num2cell(1:size(found,1));
        found = num2cell(sortrows(cell2mat(found), sortOrd)); %sort by rows
        found(:,1) = fnames(cell2mat(found(:,1))); %replace names
        CombinedRes(realI+1:realI+size(found,1),2:end-1) = found;
        CombinedRes(realI+1:realI+size(found,1),1) = {i};
        CombinedRes(realI+1:realI+size(found,1),end) = {0};
        realI = realI+size(found,1)+1;
    end
    %currently sorts by max then sum of spectra

    toc;
    % If this search is for the Jackson Lab, use the following segment:

    %     for j = 1:length(craptxt(:,1))
    %         if isempty(strcmp(FileOutCombined{i,1},craptxt{j,1})) == 0
    %             reversefinder(i,1) = 1;
    %         end;
    %     end;
    % end
    toc;

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
    HeaderFileString(1:2) = {'Rank_Number'; 'Protein_Name'};
    HeaderFileString{end-1} = 'Containment';
    HeaderFileString{end} = 'Row_Type';

    for i = 1 : NumFilesRead

        CleanedFileName = strrep(TempStruct(i).NameForFile, '.raw_20', '_');
        CleanedFileName = strrep(CleanedFileName, '_Byonic','');
        CleanedFileName = strrep(CleanedFileName, '_Elite','');
        CleanedFileName = strrep(CleanedFileName, '_control',''); %Remove after I'm done with Janos' dataset
        CleanedFileName = strrep(CleanedFileName, '25fmol_ul_6x5spike_',''); %Remove after I'm done with Janos' dataset
        %label by dataset
        for j = 1:length(UniqueColumns)
            tval = [TempStruct(1).dat.Properties.VariableNames{UniqueColumns(j)} '_' CleanedFileName];
            HeaderFileString{2+(j-1)*NumFilesRead+totUniqueFuncs(j)-length(UniqueCombineFunctions{j})+i} = matlab.lang.makeValidName(tval);
        end
    end

    %label single columns
    for j = 1:length(SingleColumns)
        HeaderFileString{2+length(UniqueColumns)*NumFilesRead+totUniqueFuncs(end)+j} = TempStruct(1).dat.Properties.VariableNames{SingleColumns(j)};
    end

    %label max and sum
    for j = 1:length(UniqueColumns)
        for jj = 1:length(UniqueCombineFunctions{j})
            HeaderFileString{2+j*NumFilesRead+totUniqueFuncs(j)-length(UniqueCombineFunctions{j})+jj} = [func2str(UniqueCombineFunctions{j}{jj}) '_' TempStruct(1).dat.Properties.VariableNames{UniqueColumns(j)}];
        end
    end

    %Now we make the final output file for printing to a csv. we first make a
    %table rearrange it appropriately, and then have it fill in the rank number
    %values before printing_ it out.

    FinalFileOut = cell2table(FileOutFinalDeComma);
    FinalFileOut.Properties.VariableNames = HeaderFileString;

    writetable(FinalFileOut,[filename '.csv']);
    fid = fopen([filename '.json'], 'w');
    fprintf(fid, jsonencode(FinalFileOut,'ConvertInfAndNaN', false));
    fclose(fid);
%     output = table2struct(FinalFileOut);
%     save('output.mat', 'output', '-v7.3');
    toc;
end