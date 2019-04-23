function [bhv, bhv_fields] = fn_read_behav_csv(csv_fname)
% Loads target time behavioral data from given csv filename
% INPUTS:
%   csv_fname [str] - full file path to behavioral csv
% OUTPUTS:
%   bhv [struct] - data coverted to struct with no spaces in variabel name
%   bhv_fields [cell array] - names of the fields of the bhv struct

% Read csv file
fprintf('\tReading behavioral csv file: %s\n',csv_fname);
bhv_file = fopen(csv_fname, 'r');
% Fields: Total_Trial, Block, Condition, Hit, RT, Timestamp, Tolerance, Trial, Score, ITI, ITI type
bhv_fields = textscan(bhv_file, '%s %s %s %s %s %s %s %s %s %s %s', 1, 'Delimiter', ',');
bhv_data = textscan(bhv_file, '%d %d %s %d %f %f %f %d %d %f %f',...
    'Delimiter', ',', 'HeaderLines', 1);
fclose(bhv_file);

% Convert field names with spaces to underscore, add to struct
for f_ix = 1:numel(bhv_fields)
    bhv_fields{f_ix} = strrep(bhv_fields{f_ix}{1},' ','_');
    bhv.(bhv_fields{f_ix}) = bhv_data{f_ix};
end

end