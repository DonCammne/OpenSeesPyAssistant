classdef LoadData
    % Class that manages multiple functions

    
    properties(Constant)
        % delimiter used in imported info from python models
        delimiter = "##############################";
        
        % current directory and indices to manipulate it
        current_dir = pwd;
        idcs   = strfind(LoadData.current_dir, '/');
    end
    
    
    methods(Static = true)
        
        function data = ImportDataFromOutside(changed_path_and_filename, column_chosen, delimiterIn, headerlinesIn)
            % Function that use "importdata" in a personalised way for
            % importing from outside of the current directory
            
            % filename_and_changed_path : vector of str, size 2
            %   Entry 1: Changed path form the current one. One backslash at the
            %   start and end is mandatory
            %   Entry 2: The names of the .txt files to be loaded from outside
            % column_chosen : vector of int
            %   Which column(s) to take from the selected data (default: 1)
            % delimiterIn : str
            %   The delimiter to separate the lines of data (default: ' ',
            %   space)
            % headerlinesIn : int
            %   The number of lines skipped (default: 0)
            
            changed_path = char(changed_path_and_filename(1));
            filename = changed_path_and_filename(2);
            
            if nargin < 2 
                column_chosen = 1;
            end
            if nargin < 3
                delimiterIn = ' ';
            end
            if nargin < 4
                headerlinesIn = 0;
            end
            
            % Check
            if changed_path(1) ~= '/' || changed_path(end) ~= '/'
                error("There should be a backslash at the start and end of "+changed_path);
            end
            
            % Parent folder
            if length(changed_path) == 1
                nr_delim = 2;
            else
                nr_delim = length(strfind(changed_path, '/'));
            end
            parent_dir = LoadData.current_dir(1:LoadData.idcs(end)-nr_delim+1);
            
            % Complete path
            path_new = strcat(parent_dir, changed_path, filename);
            
            % Import Data
            tmp = importdata(path_new, delimiterIn, headerlinesIn);
            data = tmp(:, column_chosen);
        end
        
        function data = Results(filename)
            % Function that loads the results using a list of string w/
            % the filenames in filename
            
            % filename : list of str
            %   The names of the .txt files to be loaded
            
            disp("Loading the results...");
            data = {};
            for i = 1:size(filename, 1)
                data = [data readmatrix(filename(i))];
            end
            disp("Results loaded successfully");
            disp(newline);
        end
        
        function result = ChooseResult(data, data_chosen, column_chosen)
            % Function that return a matrix with the column(s) chosen
            
            % data : cell
            %   Loaded data
            % data_chosen : int
            %   Which cell to take
            % column_chosen : vector of int
            %   Which column(s) to take from the selected cell (default: 1)
            
            if nargin < 3
                column_chosen = 1;
            end
            
            tmp = cell2mat(data(data_chosen));
            result = [0; tmp(:, column_chosen)];
        end
        
        % New
        
        function infos = InitilizeLoading(filename_saved_info)
            % Function that initialize the loading of infos form sections,
            % material models, fibers and member
            
            % filename_saved_info : str
            %   Name of the file .txt with the infos
            
            disp("Loading the infos...");
            f = fopen(filename_saved_info);
            infos = textscan(f, '%s %s', 'delimiter', '\t');
            fclose(f);
            disp("Infos loaded successfully");
            disp(newline);
        end
        
        % TODO: Implement in the fuction LoadNextInfo the argument of which
        % type of info you want to load, skipping the others (whitelist)
        % and the opposite (blacklist)
        
        function [data, next_line] = LoadNextInfo(infos, starting_line)
            % Function that loads data from a SectionSteelIShape
            
            % infos : cell of cells
            %   Information with the name of variable of each line and the value stored in cells to be retrieved easily
            % starting_line : int
            %   The line where to start retrieving infos
            
            line = starting_line;
            while line <= length(infos{1}) && infos{2}{line} ~= LoadData2.delimiter
                readl = infos{2}{line};
                if isnan(str2double(readl))
                    if readl(1) == '['
                        if readl(1:7) == '[list(['
                            readl = erase(readl, 'list(');
                            readl = erase(readl, ')');
                            data.(infos{1}{line}) = str2num(readl);
                        else
                            data.(infos{1}{line}) = str2num(readl);
                        end
                    else
                        data.(infos{1}{line}) = readl;
                    end
                else
                    data.(infos{1}{line}) = str2double(readl);
                end
                line = line + 1;
            end
            next_line = line + 1;
        end
        
    end
    
end





