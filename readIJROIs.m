function ROIs = readIJROIs(file_name_and_path)

if endsWith(file_name_and_path,'.zip')
    fnames = unzip(file_name_and_path, '_temp_unzip');
    ROIs = ROI.empty(numel(fnames),0);
    for i = 1:numel(fnames)
        ROIs(i) = ROI(fnames{i});
    end
    rmdir('_temp_unzip','s');
elseif endsWith(file_name_and_path,'.roi')
    ROIs = ROI(file_name_and_path);
else
    error('Expected a .roi file or a .zip file containing .roi files');
end

end