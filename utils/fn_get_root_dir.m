function [root_dir, app_dir] = fn_get_root_dir()
%% Check with user/OS/root directory, return relevant locations
if exist('/home/knight/hoycw/','dir')
    root_dir = '/home/knight/hoycw/';
    app_dir   = [root_dir 'Apps/'];
elseif exist('/Volumes/hoycw_clust/','dir')
    root_dir = '/Volumes/hoycw_clust/';
    app_dir   = '/Users/colinhoy/Code/Apps/';
elseif exist('/Users/SCS22/','dir');
    root_dir='/Users/SCS22/Desktop/Knight_Lab/';
    ft_dir='/Users/SCS22/Documents/MATLAB/';
else
    error('root directory not found. where are you running this?');
end

end
