
%% EEG 09 Processing Variables
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip/';
elseif exist('/Users/aasthashah/', 'dir'); root_dir = '/Users/aasthashah/Desktop/', ft_dir = '/Users/aasthashah/Applications/fieldtrip';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end


addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath(ft_dir);
ft_defaults

%% Questionares
%{ 
6 total questions:
1) How would you feel about this feedback? (Correct Easy (green) Response) 1 =
Terrible, 9 = Great
2) How would you feel about this feedback? (Incorrect Easy (red) Response) 1 =
Terrible, 9 = Great
3) How would you feel about this feedback? (Surprise Easy (blue) Response) 1 =
Terrible, 9 = Great
4) How would you feel about this feedback? (Correct Hard (green) Response) 1 =
Terrible, 9 = Great
5) How would you feel about this feedback? (Incorrect Hard (red) Response) 1 =
Terrible, 9 = Great
6) How would you feel about this feedback? (Surprise Hard (blue) Response) 1 =
Terrible, 9 = Great
%}
question_answers = {7,3,1,7,4,1};
% Response in answer to question: does the blue feedback reflect anything
% about your performance?
free_response = {'Thought it was blue when she got it right'};