classdef LoadData
    % Class that manage multiple functions

    
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
            % importing from outside of the corrent directory
            
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
            % Function that loads the ruesults using a list of string w/
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
            % Function that initialize the loading of infos form sections
            % and material models
            
            % filename_saved_info : str
            %   Name of the file .txt with the infos
            
            disp("Loading the infos...");
            f = fopen(filename_saved_info);
            infos = textscan(f, '%s %s', 'delimiter', '\t');
            fclose(f);
            disp("Infos loaded successfully");
            disp(newline);
        end
        
        %TODO: Implement function that look at the info_name and create a
        %vector with the corresponding value for each element in the info 
        
        %TODO: Implement function that returns the number of info types
        %stored in info
        
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
        
        % Old
        
        function [tlines, next_line] = InitilizeInfoLoading(filename_saved_info)
            % Function that initialize the loading of infos form sections
            % and material models
            
            % filename_saved_info : str
            %   Name of the file .txt with the infos
            
            disp("Loading the infos...");
            fid = fopen(filename_saved_info);

            tline = fgetl(fid); % Skip first line (null)
            tline = fgetl(fid);
            tlines = cell(0, 1);
            while ischar(tline)
                tlines{end+1, 1} = tline;
                tline = fgetl(fid);
            end
            fclose(fid);
            next_line = 1;
            disp("Infos loaded successfully");
            disp(newline);
        end
        
        function [data, next_line] = SectionSteelIShape(tlines, starting_line)
            % Function that loads data from a SectionSteelIShape
            
            % tlines : cells
            %   Information stored in cells to be retrieved easily
            % starting_line : int
            %   The line where to start retrieving infos
            
            disp("Loading...");
            disp(tlines{starting_line});
            disp(tlines{starting_line+1});
            data.TAG = tlines{starting_line+1};
            data.type = tlines{starting_line+2};
            data.E = str2double(tlines{starting_line+3});
            data.Fy = str2double(tlines{starting_line+4});
            data.d = str2double(tlines{starting_line+5});
            data.bf = str2double(tlines{starting_line+6});
            data.tf = str2double(tlines{starting_line+7}); 
            data.tw = str2double(tlines{starting_line+8}); 
            data.L = str2double(tlines{starting_line+9}); 
            data.r = str2double(tlines{starting_line+10}); 
            data.h_1 = str2double(tlines{starting_line+11}); 
            data.A = str2double(tlines{starting_line+12});
            data.Npl = str2double(tlines{starting_line+13});
            data.Iy = str2double(tlines{starting_line+14});
            data.Iz = str2double(tlines{starting_line+15});
            data.Wply = str2double(tlines{starting_line+16});
            data.My = str2double(tlines{starting_line+17});
            data.Iy_mod = str2double(tlines{starting_line+18});
            data.iz = str2double(tlines{starting_line+19});
            disp("Data Loaded successfully");
            disp(newline);
            
            next_line = 19 + 2 + starting_line;
        end
        
        function [data, next_line] = IMKMaterialModel(tlines, starting_line)
            % Function that loads data from a IMKMaterialModel
            
            % tlines : cells
            %   Information stored in cells to be retrieved easily
            % starting_line : int
            %   The line where to start retrieving infos
            
            disp("Loading...");
            disp(tlines{starting_line});
            disp("ID = " + tlines{starting_line+1});
            data.ID = str2double(tlines{starting_line+1});
            data.eleNameTAG = tlines{starting_line+2};
            data.L_b = str2double(tlines{starting_line+3});
            data.N_G = str2double(tlines{starting_line+4});
            data.L_o = str2double(tlines{starting_line+5});
            data.K_factor = str2double(tlines{starting_line+6});
            data.Ke = str2double(tlines{starting_line+7});
            data.a_s = str2double(tlines{starting_line+8});
            data.My_star = str2double(tlines{starting_line+9});
            data.Mc = str2double(tlines{starting_line+10});
            data.K = str2double(tlines{starting_line+11});
            data.theta_p = str2double(tlines{starting_line+12});
            data.theta_pc = str2double(tlines{starting_line+13});
            data.theta_u = str2double(tlines{starting_line+14});
            data.McMy = str2double(tlines{starting_line+15});
            data.rateDet = str2double(tlines{starting_line+16});
            data.a_s = str2double(tlines{starting_line+17});
            data.a = str2double(tlines{starting_line+18});
            disp("Data Loaded successfully");
            disp(newline);
            
            next_line = 18 + 2 + starting_line;
        end
        
        function [data, next_line] = SectionRCRectShape(tlines, starting_line)
            % Function that loads data from a SectionRCRectShape
            
            % tlines : cells
            %   Information stored in cells to be retrieved easily
            % starting_line : int
            %   The line where to start retrieving infos
            
            disp("Loading...");
            disp(tlines{starting_line});
            disp(tlines{starting_line+1});
            data.TAG = tlines{starting_line+1};
            data.fc = str2double(tlines{starting_line+2});
            data.Ec = str2double(tlines{starting_line+3});
            data.b = str2double(tlines{starting_line+4});
            data.bc = str2double(tlines{starting_line+5});
            data.d = str2double(tlines{starting_line+6});
            data.dc = str2double(tlines{starting_line+7});
            data.L = str2double(tlines{starting_line+8});
            data.e = str2double(tlines{starting_line+9});
            data.A = str2double(tlines{starting_line+10});
            data.Ac = str2double(tlines{starting_line+11});
            data.fy = str2double(tlines{starting_line+12});
            data.Ey = str2double(tlines{starting_line+13});
            data.wx = str2num(tlines{starting_line+14});
            data.wy = str2num(tlines{starting_line+15});
            data.D_bars = str2double(tlines{starting_line+16});
            data.nr_bars = str2double(tlines{starting_line+17});
            data.cl_bars = str2double(tlines{starting_line+18});
            data.rho_bars = str2double(tlines{starting_line+19});
            data.Ay = str2double(tlines{starting_line+20});
            data.fs = str2double(tlines{starting_line+21});
            data.Es = str2double(tlines{starting_line+22});
            data.D_hoops = str2double(tlines{starting_line+23});
            data.nr_hoops = str2double(tlines{starting_line+24});
            data.cl_hoops = str2double(tlines{starting_line+25});
            data.rho_s_x = str2double(tlines{starting_line+26});
            data.rho_s_y = str2double(tlines{starting_line+27});
            data.s = str2double(tlines{starting_line+28});
            data.Iy = str2double(tlines{starting_line+29});
            data.Iz = str2double(tlines{starting_line+30});
            disp("Data Loaded successfully");
            disp(newline);
            
            next_line = 30 + 2 + starting_line;
            
        end
        
        function [data, next_line] = Concrete04MaterialModel(tlines, starting_line)
            % Function that loads data from a Concrete04MaterialModel
            
            % tlines : cells
            %   Information stored in cells to be retrieved easily
            % starting_line : int
            %   The line where to start retrieving infos
            
            disp("Loading...");
            disp(tlines{starting_line});
            disp("ConfinedID = " + tlines{starting_line+1} + " UnconfinedID = " + tlines{starting_line+2});
            data.ConfID = str2double(tlines{starting_line+1});
            data.UnconfID = str2double(tlines{starting_line+2});
            data.eleTAG = tlines{starting_line+3};
            data.ec = str2double(tlines{starting_line+4});
            data.ecp = str2double(tlines{starting_line+5});
            data.esu = str2double(tlines{starting_line+6});
            data.ecu = str2double(tlines{starting_line+7});
            data.k1 = str2double(tlines{starting_line+8});
            data.k2 = str2double(tlines{starting_line+9});
            data.ineffectual_area = str2double(tlines{starting_line+10});
            data.Ae = str2double(tlines{starting_line+11});
            data.rho_cc = str2double(tlines{starting_line+12});
            data.Acc = str2double(tlines{starting_line+13});
            data.ke = str2double(tlines{starting_line+14});
            data.fl_x = str2double(tlines{starting_line+15});
            data.fl_y = str2double(tlines{starting_line+16});
            data.fl_prime = str2double(tlines{starting_line+17});
            data.confinement_factor = str2double(tlines{starting_line+18});
            data.fcc = str2double(tlines{starting_line+19});
            data.ecc = str2double(tlines{starting_line+20});
            data.eccu = str2double(tlines{starting_line+21});
            data.N_G = str2double(tlines{starting_line+22});
            data.L_o = str2double(tlines{starting_line+23});
            data.K_factor = str2double(tlines{starting_line+24});
            disp("Data Loaded successfully");
            disp(newline);
            
            next_line = 24 + 2 + starting_line;
            
        end
        
        function [data, next_line] = Steel01forRCMaterialModel(tlines, starting_line)
            % Function that loads data from a Steel01forRCMaterialModel
            
            % tlines : cells
            %   Information stored in cells to be retrieved easily
            % starting_line : int
            %   The line where to start retrieving infos
            
            disp("Loading...");
            disp(tlines{starting_line});
            disp("ID = " + tlines{starting_line+1});
            data.ID = str2double(tlines{starting_line+1}); 
            data.eleNameTAG = tlines{starting_line+2}; 
            data.b = str2double(tlines{starting_line+3}); 
            disp("Data Loaded successfully");
            disp(newline);
            
            next_line = 3 + 2 + starting_line;
            
        end
        
        function [data, next_line] = FibersRectRCShape(tlines, starting_line)
            % Function that loads data from a FibersRectRCShape
            
            % tlines : cells
            %   Information stored in cells to be retrieved easily
            % starting_line : int
            %   The line where to start retrieving infos
            
            disp("Loading...");
            disp(tlines{starting_line});
            disp("ID = " + tlines{starting_line+1});
            data.ID = str2double(tlines{starting_line+1});
            data.eleNameTAG = tlines{starting_line+2};
            data.confID = str2double(tlines{starting_line+3});
            data.unconfID = str2double(tlines{starting_line+4});
            data.barsID = str2double(tlines{starting_line+5});
            data.discr_core = str2num(tlines{starting_line+6});
            data.discr_cover_lateral = str2num(tlines{starting_line+7});
            data.discr_cover_updown = str2num(tlines{starting_line+8});
            data.xpos_history = str2num(tlines{starting_line+9});
            data.ypos_history = str2num(tlines{starting_line+10});
            disp("Data Loaded successfully");
            disp(newline);
            
            next_line = 10 + 2 + starting_line;
            
        end
        
        function [data, next_line] = PZRotSpringMaterialModel(tlines, starting_line)
            % Function that loads data from a PZRotSpringMaterialModel
            
            % tlines : cells
            %   Information stored in cells to be retrieved easily
            % starting_line : int
            %   The line where to start retrieving infos
            
            disp("Loading...");
            disp(tlines{starting_line});
            disp("ID  = " + tlines{starting_line+1});
            data.ConfID = str2double(tlines{starting_line+1});
            data.colNameTAG = tlines{starting_line+2};
            data.d_beam = str2double(tlines{starting_line+3});
            data.a_s = str2double(tlines{starting_line+4});
            data.Ry = str2double(tlines{starting_line+5});
            data.pinchx = str2double(tlines{starting_line+6});
            data.pinchy = str2double(tlines{starting_line+7});
            data.dmg1 = str2double(tlines{starting_line+8});
            data.dmg2 = str2double(tlines{starting_line+9});
            data.beta = str2double(tlines{starting_line+10});
            data.Vy = str2double(tlines{starting_line+11});
            data.G = str2double(tlines{starting_line+12});
            data.Ke = str2double(tlines{starting_line+13});
            data.Kp = str2double(tlines{starting_line+14});
            data.gamma1_y = str2double(tlines{starting_line+15});
            data.gamma2_y = str2double(tlines{starting_line+16});
            data.gamma3_y = str2double(tlines{starting_line+17});
            data.M1y = str2double(tlines{starting_line+18});
            data.M2y = str2double(tlines{starting_line+19});
            data.M3y = str2double(tlines{starting_line+20});
            data.Ks_ref = str2double(tlines{starting_line+20});   % refined parameters
            data.Kb_ref = str2double(tlines{starting_line+21});
            data.Ke_ref = str2double(tlines{starting_line+22});
            data.Ksf = str2double(tlines{starting_line+23});
            data.Kbf = str2double(tlines{starting_line+24});
            data.Kf = str2double(tlines{starting_line+25});
            data.Kf_Ke = str2double(tlines{starting_line+26});
            data.Cw1 = str2double(tlines{starting_line+27});
            data.Cf1 = str2double(tlines{starting_line+28});
            data.Cw4 = str2double(tlines{starting_line+29});
            data.Cf4 = str2double(tlines{starting_line+30});
            data.Cw6 = str2double(tlines{starting_line+31});
            data.Cf6 = str2double(tlines{starting_line+32});
            data.V1 = str2double(tlines{starting_line+33});
            data.V4 = str2double(tlines{starting_line+34});
            data.V6 = str2double(tlines{starting_line+35});
            data.M1 = str2double(tlines{starting_line+36});
            data.M4 = str2double(tlines{starting_line+37});
            data.M6 = str2double(tlines{starting_line+38});
            data.Gamma_1 = str2double(tlines{starting_line+39});
            data.Gamma_4 = str2double(tlines{starting_line+40});
            data.Gamma_6 = str2double(tlines{starting_line+41});
            disp("Data Loaded successfully");
            disp(newline);
            
            next_line = 41 + 2 + starting_line;
            
        end
        
    end
    
    
end





