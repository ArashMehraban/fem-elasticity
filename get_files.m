function files = get_files(folderName, file_name_list)
% function files = get_files(all_files, file_name_list)
%
%GET_FILES returns mesh files stored in folderName folder according to
%the list of file names specified in file_name_list
%
%input:     folderName: the folder where mesh files are stored
%     : file_name_list: list of filenames to be extracted from folderName
%                        
%output: mesh files from folderName specified in file_name_list



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