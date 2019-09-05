function [root_dir, app_dir] = fn_get_root_dir()
%% Check with user/OS/root directory, return relevant locations
if exist('/home/knight/','dir')
    root_dir = '/home/knight/';
    app_dir   = [root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Volumes/hoycw_clust/','dir')
    root_dir = '/Volumes/hoycw_clust/';
    app_dir   = '/Users/colinhoy/Code/Apps/';
elseif exist('/Users/sheilasteiner/','dir');
    root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';
    app_dir='/Users/sheilasteiner/Documents/MATLAB/';
else
    error('root directory not found. where are you running this?');
end

end
