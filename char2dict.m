function Dictionary = char2dict(charArray)
cellArray = strsplit(charArray, newline);

Dictionary = containers.Map();

for i = 1:length(cellArray)
    if strfind(cellArray{i}, ' = ')
        splitCells = strsplit(cellArray{i}, ' = ');
        num = str2double(splitCells{2});
        
        if ~isnan(num)
            splitCells{2} = num;
        end
        Dictionary(splitCells{1}) = splitCells{2};
    end
end
end