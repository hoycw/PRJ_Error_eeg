
%% EEG Pilot 08 Processing Variables
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip/';
elseif exist('/Users/aasthashah/', 'dir'); root_dir = '/Users/aasthashah/Desktop/', ft_dir = '/Users/aasthashah/Applications/fieldtrip';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end


addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath(ft_dir);
ft_defaults

%% Questionares
%{
    12 total questions:
    1) How would you feel about this feedback? (Correct Easy (green) Response) 1 =
    Terrible, 9 = Great
    2) Does this feedback tell you about your performance 1 = Yes 0 = No 2 = I
    don't know
    3) How would you feel about this feedback? (Incorrect Easy (red) Response) 1 =
    Terrible, 9 = Great
    4) Does this feedback tell you about your performance 1 = Yes 0 = No 2 = I
    don't know
    5) How would you feel about this feedback? (Surprise Easy (blue) Response) 1 =
    Terrible, 9 = Great
    6) Does this feedback tell you about your performance 1 = Yes 0 = No 2 = I
    don't know
    7) How would you feel about this feedback? (Correct Hard (green) Response) 1 =
    Terrible, 9 = Great
    8) Does this feedback tell you about your performance 1 = Yes 0 = No 2 = I
    don't know
    9) How would you feel about this feedback? (Incorrect Hard (red) Response) 1 =
    Terrible, 9 = Great
    10) Does this feedback tell you about your performance 1 = Yes 0 = No 2 = I
    don't know
    11) How would you feel about this feedback? (Surprise Hard (blue) Response) 1 =
    Terrible, 9 = Great
    12) Does this feedback tell you about your performance 1 = Yes 0 = No 2 = I
    don't know
    %}
    question_answers = {9, 2, 5, 9, 6, 5};
    % Response in answer to qustion: does the blue feedback reflect anything
    % about your performance?
    free_response = {'Thought it was random'};
    
    

