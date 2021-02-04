function [aud_sal] = fn_load_auditory_salience(csv_fname)
%% Load auditory salience info per sound from csv
%   Auditory salience metrics per sound provided by Sijia Zhao
% INPUTS:
%   csv_fname [str] - full file path to behavioral csv
% OUTPUTS:
%   aud_sal [struct] - contains 6 metrics of auditory salience per sound
% CSV format:
%   filename [str] - name of the sound file (e.g., "P3A01.WAV")
%   ERB [float] - "equivalent rectangular bandwidth", roughly loudness
%       "Loudness" is the strongest driver of salience!
%   roughness [float] - see Zhao 2019 JNeuro; strongly correlates with subjective
%       salience and physiological index of saleince (pupil dilation)
%   maxSalience [float] - max amplitude in Kayser (2006) salience map model
%   meanSalience [float] - mean amplitude in Kayser (2006) salience map model
%   maxGradient [float] -  max change in time in Kayser (2006) salience map model
%   meanGradient [float] - mean change in time in Kayser (2006) salience map model


% Open file
fprintf('\tReading behavioral csv file: %s\n',csv_fname);
bhv_file = fopen(csv_fname);

% Read field names
csv_fields  = {'filename', 'ERB', 'roughness', 'maxSalience', 'meanSalience', 'maxGradient', 'meanGradient'};
out_fields = {'sound', 'ERB', 'rough', 'mxS', 'mnS', 'mxdS', 'mndS'};
file_cols = textscan(bhv_file, '%s %s %s %s %s %s %s', 1, 'Delimiter', ',');

% Check that loaded fields match expected py_fields
if numel(csv_fields)~=numel(file_cols)
    error('Mismatched fields in auditory salience csv and expected!');
end
for f = 1:numel(csv_fields)
   if ~strcmp(csv_fields{f},file_cols{f}{1})
      error(['Mismatched field in behav csv and expected: ' csv_fields{f} ' vs. ' file_cols{f}{1}]);
    end
end

% Read data
sal_data = textscan(bhv_file,'%s %f %f %f %f %f %f','Delimiter',',','HeaderLines',1);
fclose(bhv_file);
n_sounds = size(sal_data{1},1);
fprintf('\t\tFound %d sounds in auditory salience file\n', n_sounds);

% Add data to struct
for f_ix = 1:numel(file_cols)
    if numel(sal_data{f_ix})==n_sounds
        aud_sal.(out_fields{f_ix}) = sal_data{f_ix};
    else
        error([csv_fname ' column ' file_cols{f_ix} ' has unexpected number of rows']);
    end
end

end
