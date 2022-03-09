function outputStruct = parseScanImageMetaData(inputStr)
    inputStr = string(inputStr);
    inputStr = splitlines(inputStr);

    fieldNames = [];
    values = {};
    
    allNames = extractBefore(inputStr, " = ");
    splitInds = regexp(allNames, "\D\.\D", 'once', 'forceCellOutput');
    needSplit = ~cellfun(@isempty, splitInds);
    
    % save all the lines that do not need splitting
    readyForValStr = inputStr(~needSplit);
    for i = 1:length(readyForValStr)
        nameAndVal = split(readyForValStr(i), " = ");
        if length(nameAndVal) == 2
            fieldNames = [fieldNames; nameAndVal(1)];
            
            val = nameAndVal(2);
            number = str2num(val);

            if isempty(number) || isa(number,'function_handle')
                values = [values; {val}];
            else
                values = [values; number];
            end      
        end
    end
    
    % find the ones that need further splitting and recursively return
    if any(needSplit)
        needSplitStr = inputStr(needSplit);
        splitIndsAarr = [splitInds{needSplit}]' + 1;
        structNames = extractBefore(needSplitStr, splitIndsAarr);
        remainingString = extractAfter(needSplitStr, splitIndsAarr);   
        uniqueNames = unique(structNames);
        for i = 1:length(uniqueNames)
            newStringForStruct = remainingString(structNames == uniqueNames(i));
            subStruct = parseScanImageMetaData(newStringForStruct);

            fieldNames = [fieldNames; uniqueNames(i)];
            values = [values; subStruct];
        end
    end

    outputStruct = cell2struct(values, fieldNames);
end

