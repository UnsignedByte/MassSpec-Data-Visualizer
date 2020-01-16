function getPeptideMap(proteinName, dat, summarydat, fastaFile, filename, sheetname)
    modlist = getAllMods(summarydat);
    Istart = 0;

    for j = 1:size(dat, 1) %get range of rows with specified protein rank
        if strcmp(dat.ProteinName(j), proteinName) && Istart == 0
            Istart = j;
        elseif ~strcmp(dat.ProteinName(j), proteinName) && Istart ~= 0
            Iend = j-1;
            break;
        end
    end

    proteinSequence = NaN;
    
    for i = 1:size(fastaFile,1) %find protein sequence from fasta
        if strcmp(fastaFile.Header{i},proteinName(2:end)) %2:end used to ignore the preceding >
            proteinSequence = fastaFile.Sequence{i};
            break;
        end
    end
    
    if isnan(proteinSequence)
        disp(['protein ' proteinName ' not found']);
        return;
    end

    tablewidth = size(modlist, 1)+3;

    alltypes    = cell(tablewidth, 1); %types of the thing
    alltypes(:) = {'uint32'};
    alltypes{2} = 'char';
    allvarnames = cell(tablewidth, 1);
    allvarnames(1:3) = {'Position', 'AA', 'Total_Peptides'};
    allvarnames(4:end) = modlist(:,1);

    resTable = table('Size', [length(proteinSequence), tablewidth],'VariableTypes',alltypes, 'VariableNames', allvarnames);

    dqueue = zeros((Iend-Istart+1)*2,3);
    peps = cell(2,Iend-Istart+1);

    for i = Istart:Iend
        [pep, mods] = formatPeptide(dat.Peptide_ProteinMetricsConfidential_{i}, modlist);
        peps{1,i-Istart+1} = pep;
        peps{2,i-Istart+1} = mods;
        dqueue((i-Istart)*2+(1:2),:) = [dat.StartingPosition(i), i-Istart+1, 1; dat.StartingPosition(i)+length(pep)-4, i-Istart+1, -1];
    end
    dqueue = sortrows(dqueue);

    for i = 1:length(proteinSequence)
        if size(dqueue,1) > 0
            lastSlice = resTable(max(i-1,1),3);
            while dqueue(1,1) == i
                lastSlice.Total_Peptides(1) = lastSlice.Total_Peptides(1)+dqueue(1,3);
                for j = 1:size(peps{2,dqueue(1,2)},1)
        %             disp(peps{2,dqueue(1,2)});
                    if dqueue(1,3) == 1
                        resTable.(peps{2,dqueue(1,2)}{j,1})(max(i+peps{2, dqueue(1,2)}{j,2}-3,1)) = resTable.(peps{2,dqueue(1,2)}{j,1})(max(i+peps{2, dqueue(1,2)}{j,2}-3,1)) + 1;
                    end
        %             lastSlice.(peps{2,dqueue(1,2)}{j,1})(1) = lastSlice.(peps{2,dqueue(1,2)}{j,1})(1)+dqueue(1,3);
                end
                dqueue(1,:) = [];
                if size(dqueue,1) == 0
                    break;
                end
            end
            resTable(i,3) = lastSlice;
        end
    %     disp(lastSlice);
        resTable.Position(i) = i;
        resTable.AA{i} = proteinSequence(i);
    end
    todel = [];
    for i = 4:size(resTable,2)
        if ~any(table2array(resTable(:,i)))
            todel(end+1) = i;
        end
    end
    resTable(:,todel) = [];

    writetable(resTable,[filename '.xlsx'], 'Sheet', sheetname);
end