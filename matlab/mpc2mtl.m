function v = mpc2mtl(filename)
fid = fopen(filename, 'r');

%v.b = zeros(16,10000);

a = 0;
%countSubj = 1;
while ~strcmp(a,'')
    a = fscanf(fid, '%s',1);
    
    % --- the following part will be used in the future ---
    if strcmp(a,'Subject:')        
        aux  = fscanf(fid,'%s',1);                  % scans subject name/number 
        subj = str2double(aux);                     % converts to number and stores in subj
        if isnan(subj)                              % checks if the label was converted to number
            nPref = find_prefix(aux);
            fprintf(1,'%s',aux(1:nPref-1));         % print the prefix
            subj = str2double(aux(nPref:end));      % removes the prefix
        end        
        if subj == 0 || isnan(subj)
            warning('subject = 0 or NaN found, it will be skipped \n');
        else
            fprintf(1,'%d,  ',subj);                    % print on screen
            v.s(subj) = subj;                           % stores in field s
        end
    end

    if exist('subj','var') && subj ~= 0 && ~isnan(subj)
        % --- reading var B ---
        if strcmp(a,'B:')
            for i=1:5:10000-5
                a = fscanf(fid, '%s',1);
                v.b(i  ,subj) = str2double(fscanf(fid, '%s',1));
                v.b(i+1,subj) = str2double(fscanf(fid, '%s',1));
                v.b(i+2,subj) = str2double(fscanf(fid, '%s',1));
                v.b(i+3,subj) = str2double(fscanf(fid, '%s',1));
                v.b(i+4,subj) = str2double(fscanf(fid, '%s',1));
            end
        end

        % --- reading var G ---
        if strcmp(a,'G:')
            for i=1:5:1000-5
                a = fscanf(fid, '%s',1);
                v.g(i  ,subj) = str2double(fscanf(fid, '%s',1));
                v.g(i+1,subj) = str2double(fscanf(fid, '%s',1));
                v.g(i+2,subj) = str2double(fscanf(fid, '%s',1));
                v.g(i+3,subj) = str2double(fscanf(fid, '%s',1));
                v.g(i+4,subj) = str2double(fscanf(fid, '%s',1));
            end
        end

        % --- reading var Q ---
        if strcmp(a,'Q:')
            for i=1:5:5000-5
                a = fscanf(fid, '%s',1);
                v.q(i  ,subj) = str2double(fscanf(fid, '%s',1));
                v.q(i+1,subj) = str2double(fscanf(fid, '%s',1));
                v.q(i+2,subj) = str2double(fscanf(fid, '%s',1));
                v.q(i+3,subj) = str2double(fscanf(fid, '%s',1));
                v.q(i+4,subj) = str2double(fscanf(fid, '%s',1));
            end
        end
    end
end
fclose(fid);

%% function that looks for a prefix in the subject name
function n = find_prefix(name)
n = 2;
while isnan(str2double(name(n:end))) && n<=length(name)
    n = n+1;
end


    





