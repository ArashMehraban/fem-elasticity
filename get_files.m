function files = get_files(folderName, file_name_list)

    files_in_folder = dir(folderName);
    all_files = {files_in_folder.name};
    
    files = cell(1,size(file_name_list,2));
    k=1;
    for i=1:size(file_name_list,2)
        files{k} = validatestring(strcat(file_name_list{i},'.exo'),all_files);
        k=k+1;
     end
    
    for i=1:size(files,2)
        files{i}= files{i}(1:end-4);
    end

end