% %EFF21 is comment signaling change to make code more efficient
% use single for spike_val but not for spike_times
% spikes_trigger and spikes_trigger_stop are not yet sparse
% update sz_on_index and sz_off_index in own function and check for
% emptyness before trying to save
% spikes_trigger to sparse? spikes_trigger initialization
% when loading data and then resaving them it might be slightly off with
% being exactly every 600 time (60*10), so saving vars won't quite work.
% when to save if only running part of file but dont push stop
% spike_count_history becomes [spike_count_history spike_count_raw
% spike_count_combined] to reduce number of times we save to file. (columns
% next to each other.

%% INITIALIZATION OF GUI
function varargout = Bayes_Op_estim_2021_09_01(varargin)
% Bayes_Op_estim_2021_09_01 MATLAB code for Bayes_Op_estim_2021_09_01.fig
%      Bayes_Op_estim_2021_09_01, by itself, creates a new Bayes_Op_estim_2021_09_01 or raises the existing
%      singleton*.
%
%      H = Bayes_Op_estim_2021_09_01 returns the handle to a new Bayes_Op_estim_2021_09_01 or the handle to
%      the existing singleton*.
%
%      Bayes_Op_estim_2021_09_01('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in Bayes_Op_estim_2021_09_01.M with the given input arguments.
%
%      Bayes_Op_estim_2021_09_01('Property','Value',...) creates a new Bayes_Op_estim_2021_09_01 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Bayes_Op_estim_2021_09_01_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Bayes_Op_estim_2021_09_01_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Bayes_Op_estim_2021_09_01

% Last Modified by GUIDE v2.5 01-Sep-2021 22:26:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Bayes_Op_estim_2021_09_01_OpeningFcn, ...
                   'gui_OutputFcn',  @Bayes_Op_estim_2021_09_01_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

%% --- Executes just before Bayes_Op_estim_2021_09_01 is made visible.
function Bayes_Op_estim_2021_09_01_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Bayes_Op_estim_2021_09_01 (see VARARGIN)

% Choose default command line output for bayes_opt_estim_v62

%Clear command window
clc 

%Reset the Data Aquisition and Image Aquisition
daqreset
if license('test','image_acquisition_toolbox')
    try
        imaqreset
    catch
        disp('Cannot run imaqreset - computer might not have image aquisition toolbox'); 
        warndlg('Cannot run imaqreset - computer might not have image aquisition toolbox'); 
    end
end

handles.output = hObject;

%initialize flags used in offline detection and in troubleshooting
global DoTroubleshoot DoRestartOffline DoStopOffline DoPauseOffline
global DoLog
global DoSaveToFile DoResaveDataFiles DoFinishedFile
DoTroubleshoot=0;
DoStopOffline=0;
DoRestartOffline=0;
DoPauseOffline=0;
DoResaveDataFiles=0;    %in offline mode, this flag will force resaving of individual data files ???!!!???built in as option?
DoLog=2;                %1=create 1 log file for entire recording, 2=create log for every hour of running,0=no log file
DoSaveToFile=1;         %standard is that program saves output to files (might get overridden by offline detection
DoFinishedFile=0;      %for offline reanalysis, we have to know when to start new file
% Update handles structured
guidata(hObject, handles);

%% --- Outputs from this function are returned to the command line.
function varargout = Bayes_Op_estim_2021_09_01_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% FUNCTIONS CALLED WHEN PRESSING A BUTTON
% =====================================================================

%% --- Executes on button press: LOAD SETTINGS
% Will load the settings file provided in the input box next to it
function PB_LoadSettings_Callback(hObject, eventdata, handles)
% hObject    handle to PB_LoadSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global DoTroubleshoot
if DoTroubleshoot
   disp('Run: PB_LoadSettings_Callback...')
end
%Load in the settings from the file provided
load(get(handles.ET_load_path,'String'));

%Overwrite the values from the structure struc to the GUI
% try
    set(handles.ET_spike_dist_min,'String',struc.ET_spike_dist_min);
% end
% try
    set(handles.ET_spike_dist_max,'String',struc.ET_spike_dist_max);
% end
% try
    set(handles.ET_width_at_percent_height,'String',struc.ET_width_at_percent_height);
% end
% try
    set(handles.ET_threshold_pos,'String',struc.ET_threshold_pos);
% end
% try
    set(handles.ET_threshold_neg,'String',struc.ET_threshold_neg);
% end
% try
    set(handles.ET_maxAmp_pos,'String',struc.ET_maxAmp_pos);
% end
% try
    set(handles.ET_maxAmp_neg,'String',struc.ET_maxAmp_neg);
% end
% try
    set(handles.ET_min_width,'String',struc.ET_min_width);
% end
% try
    set(handles.ET_max_width,'String',struc.ET_max_width);
% end
% try
    set(handles.ET_spike_dist_min_raw,'String',struc.ET_spike_dist_min_raw);
% end
% try
    set(handles.ET_spike_dist_max_raw,'String',struc.ET_spike_dist_max_raw);
% end
% try
    set(handles.ET_width_at_percent_height_raw,'String',struc.ET_width_at_percent_height_raw);
% end
% try
    set(handles.ET_threshold_pos_raw,'String',struc.ET_threshold_pos_raw);
% end
% try
    set(handles.ET_threshold_neg_raw,'String',struc.ET_threshold_neg_raw);
% end
% try
    set(handles.ET_maxAmp_raw_pos,'String',struc.ET_maxAmp_raw_pos);
% end
% try
    set(handles.ET_maxAmp_raw_neg,'String',struc.ET_maxAmp_raw_neg);
% end
% try
    set(handles.ET_min_width_raw,'String',struc.ET_min_width_raw);
% end
% try
    set(handles.ET_max_width_raw,'String',struc.ET_max_width_raw);
% end
% try
    set(handles.ET_min_spikes_per_two_s,'String',struc.ET_min_spikes_per_two_s);
% end
% try
    set(handles.ET_min_spikes_per_two_s_stop,'String',struc.ET_min_spikes_per_two_s_stop);
% end
% try
    set(handles.ET_spike_logic_1,'String',struc.ET_spike_logic_1);
% end
% try
    set(handles.ET_spike_logic_2,'String',struc.ET_spike_logic_2);
% end
% try
    set(handles.ET_spike_logic_3,'String',struc.ET_spike_logic_3);
% end
% try
    set(handles.ET_spike_logic_4,'String',struc.ET_spike_logic_4);
% end
% try
    set(handles.ET_exclusion_time,'String',struc.ET_exclusion_time);
% end
% try
    set(handles.ET_and_window,'String',struc.ET_and_window);
% end
% try
    set(handles.ET_n_ch_in,'String',struc.ET_n_ch_in);
% end
% try
    set(handles.ET_n_ch_out,'String',struc.ET_n_ch_out);
% end
% try
    set(handles.ET_seizure_detection_channels,'String',struc.ET_seizure_detection_channels);
% end
% try
    set(handles.ET_n_cams,'String',struc.ET_n_cams);
% end
% try
    set(handles.ET_fs,'String',struc.ET_fs);
% end
% try
    set(handles.ET_MonoOrBiPhasic,'String',struc.ET_MonoOrBiPhasic);
% end
% try
    set(handles.ET_ampLo,'String',struc.ET_ampLo);
% end
% try
    set(handles.ET_widthLo,'String',struc.ET_widthLo);
% end
% try
    set(handles.ET_freqLo,'String',struc.ET_freqLo);
% end
% try
    set(handles.ET_TrainDuration,'String',struc.ET_TrainDuration);
% end
% try
    set(handles.ET_SaveFolder,'String',struc.ET_SaveFolder);
% end
% try
    set(handles.ET_SaveName,'String',struc.ET_SaveName);
% end
% try
    set(handles.ET_device_ID,'String',struc.ET_device_ID);
% end
% try
    set(handles.ET_time_plot,'String',struc.ET_time_plot);
% end
% try
    set(handles.ET_exclusion_time,'String',struc.ET_exclusion_time);
% end
    set(handles.ET_channel_spacing,'String',struc.ET_channel_spacing); 
% try
    set(handles.ET_stim_dev,'String',struc.ET_stim_dev);
% end
% try 
    set(handles.ET_fs_stim,'String',struc.ET_fs_stim);
% end
% try 
    set(handles.ET_chgRatioHi,'String',struc.ET_chgRatioHi);
% end
% try 
    set(handles.ET_ampHi,'String',struc.ET_ampHi);
% end
% try 
    set(handles.ET_freqHi,'String',struc.ET_freqHi);
% end
% try 
    set(handles.ET_widthHi,'String',struc.ET_widthHi);
% end
% try
    set(handles.ET_path_for_restart,'String',struc.ET_path_for_restart);
% end
% try
    set(handles.ET_widthBins,'String',struc.ET_ampBins)
% end
% try
    set(handles.ET_widthBins,'String',struc.ET_freqBins)
% end
% try
    set(handles.ET_widthBins,'String',struc.ET_widthBins)
% end
% try 
    set(handles.ET_ampRatioLo,'String',struc.ET_ampRatioLo);
% end
% try
    set(handles.ET_ampRatioHi,'String',struc.ET_ampRatioHi)
% end
% try
    set(handles.ET_ampRatioBins,'String',struc.ET_ampRatioBins)
% end

% try 
    set(handles.ET_chgRatioLo,'String',struc.ET_chgRatioLo);
% end
% try
    set(handles.ET_chgRatioHi,'String',struc.ET_chgRatioHi)
% end
% try
    set(handles.ET_chgRatioBins,'String',struc.ET_chgRatioBins)
% end
try
    set(handles.ET_optimize_vec,'String',struc.ET_optimize_vec)
end
try
    set(handles.ET_non_optimized_val_vec,'String',struc.ET_non_optimized_val_vec)
end
try
    set(handle.ET_non_optimized_val_vec_ch1,'String',struc.ET_non_optimized_val_vec_ch1)
end
try
    set(handle.ET_non_optimized_val_vec_ch2,'String',struc.ET_non_optimized_val_vec_ch2)
end
try
    set(handle.ET_non_optimized_val_vec_ch3,'String',struc.ET_non_optimized_val_vec_ch3)
end
try
    set(handle.ET_non_optimized_val_vec_ch4,'String',struc.ET_non_optimized_val_vec_ch4)
end
try
    set(handle.ET_non_optimized_val_vec,'String',struc.ET_non_optimized_val_vec)
end

try
    set(handle.ET_non_optimized_val_vec_ch1_worst,'String',struc.ET_non_optimized_val_vec_ch1_worst)
end
try
    set(handle.ET_non_optimized_val_vec_ch2_worst,'String',struc.ET_non_optimized_val_vec_ch2_worst)
end
try
    set(handle.ET_non_optimized_val_vec_ch3_worst,'String',struc.ET_non_optimized_val_vec_ch3_worst)
end
try
    set(handle.ET_non_optimized_val_vec_ch4_worst,'String',struc.ET_non_optimized_val_vec_ch4_worst)
end
guidata(hObject,handles);


%% --- SAVE_SETTINGS is called by PB_Go_Callback (after pressing GO)
function save_settings(path, settings_name, handles)
% hObject    handle to PB_SaveSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global DoTroubleshoot
if DoTroubleshoot
   disp('Run: save_settings...')
end

struc.ET_n_ch_in=get(handles.ET_n_ch_in,'String');
struc.ET_n_ch_out=get(handles.ET_n_ch_out,'String');
struc.ET_seizure_detection_channels=get(handles.ET_seizure_detection_channels,'String');
struc.ET_n_cams=get(handles.ET_n_cams,'String');
struc.ET_fs=get(handles.ET_fs,'String');
struc.ET_MonoOrBiPhasic=get(handles.ET_MonoOrBiPhasic,'String');
struc.ET_ampLo=get(handles.ET_ampLo,'String');
struc.ET_widthLo=get(handles.ET_widthLo,'String');
struc.ET_freqLo=get(handles.ET_freqLo,'String');
struc.ET_TrainDuration=get(handles.ET_TrainDuration,'String');
struc.ET_SaveFolder=get(handles.ET_SaveFolder,'String');
struc.ET_SaveName=get(handles.ET_SaveName,'String');
struc.ET_load_path=get(handles.ET_load_path,'String');

struc.ET_spike_dist_min=get(handles.ET_spike_dist_min,'String');
struc.ET_spike_dist_max=get(handles.ET_spike_dist_max,'String');
struc.ET_width_at_percent_height=get(handles.ET_width_at_percent_height,'String');
struc.ET_threshold_pos=get(handles.ET_threshold_pos,'String');
struc.ET_threshold_neg=get(handles.ET_threshold_neg,'String');
struc.ET_maxAmp_pos=get(handles.ET_maxAmp_pos,'String');
struc.ET_maxAmp_neg=get(handles.ET_maxAmp_neg,'String');
struc.ET_min_width=get(handles.ET_min_width,'String');
struc.ET_max_width=get(handles.ET_max_width,'String');

struc.ET_spike_dist_min_raw=get(handles.ET_spike_dist_min_raw,'String');
struc.ET_spike_dist_max_raw=get(handles.ET_spike_dist_max_raw,'String');
struc.ET_width_at_percent_height_raw=get(handles.ET_width_at_percent_height_raw,'String');
struc.ET_threshold_pos_raw=get(handles.ET_threshold_pos_raw,'String');
struc.ET_threshold_neg_raw=get(handles.ET_threshold_neg_raw,'String');
struc.ET_maxAmp_raw_pos=get(handles.ET_maxAmp_raw_pos,'String');
struc.ET_maxAmp_raw_neg=get(handles.ET_maxAmp_raw_neg,'String');
struc.ET_min_width_raw=get(handles.ET_min_width_raw,'String');
struc.ET_max_width_raw=get(handles.ET_max_width_raw,'String');

struc.ET_min_spikes_per_two_s=get(handles.ET_min_spikes_per_two_s,'String');
struc.ET_min_spikes_per_two_s_stop=get(handles.ET_min_spikes_per_two_s_stop,'String');
struc.ET_spike_logic_1=get(handles.ET_spike_logic_1,'String');
struc.ET_spike_logic_2=get(handles.ET_spike_logic_2,'String');
struc.ET_spike_logic_3=get(handles.ET_spike_logic_3,'String');
struc.ET_spike_logic_4=get(handles.ET_spike_logic_4,'String');
struc.ET_exclusion_time = get(handles.ET_exclusion_time,'String');
struc.ET_and_window = get(handles.ET_and_window,'String');

struc.ET_device_ID = get(handles.ET_device_ID,'String');
struc.ET_time_plot = get(handles.ET_time_plot,'String');
struc.ET_channel_spacing = get(handles.ET_channel_spacing,'String');
struc.ET_stim_dev = get(handles.ET_stim_dev,'String');
struc.ET_fs_stim = get(handles.ET_fs_stim,'String');
struc.ET_chgRatioHi = get(handles.ET_chgRatioHi,'String');
struc.ET_ampHi = get(handles.ET_ampHi,'String');
struc.ET_freqHi = get(handles.ET_freqHi,'String');
struc.ET_widthHi = get(handles.ET_widthHi,'String');
struc.ET_path_for_restart = get(handles.ET_path_for_restart,'String');
struc.ET_ampBins = get(handles.ET_ampBins,'String');
struc.ET_freqBins = get(handles.ET_freqBins,'String');
struc.ET_widthBins = get(handles.ET_widthBins,'String');

struc.ET_ampRatioLo =  get(handles.ET_ampRatioLo,'String');
struc.ET_ampRatioHi =  get(handles.ET_ampRatioHi,'String');
struc.ET_ampRatioBins = get(handles.ET_ampRatioBins,'String');

struc.ET_chgRatioLo =  get(handles.ET_chgRatioLo,'String');
struc.ET_chgRatioHi =  get(handles.ET_chgRatioHi,'String');
struc.ET_chgRatioBins = get(handles.ET_chgRatioBins,'String');

struc.ET_optimize_vec = get(handles.ET_optimize_vec,'String');

struc.ET_non_optimized_val_vec_ch1 = get(handles.ET_non_optimized_val_vec_ch1 ,'String');
struc.ET_non_optimized_val_vec_ch2 = get(handles.ET_non_optimized_val_vec_ch2 ,'String');
struc.ET_non_optimized_val_vec_ch3 = get(handles.ET_non_optimized_val_vec_ch3 ,'String');
struc.ET_non_optimized_val_vec_ch4 = get(handles.ET_non_optimized_val_vec_ch4 ,'String');

struc.ET_non_optimized_val_vec = get(handles.ET_non_optimized_val_vec,'String');

struc.ET_non_optimized_val_vec_ch1_worst = get(handles.ET_non_optimized_val_vec_ch1_worst ,'String');
struc.ET_non_optimized_val_vec_ch2_worst = get(handles.ET_non_optimized_val_vec_ch2_worst ,'String');
struc.ET_non_optimized_val_vec_ch3_worst = get(handles.ET_non_optimized_val_vec_ch3_worst ,'String');
struc.ET_non_optimized_val_vec_ch4_worst = get(handles.ET_non_optimized_val_vec_ch4_worst ,'String');

save([path '\' settings_name '.mat'], 'struc');


%% --- Executes on button press in PB_Go: PRESSING THE GO BUTTON
function PB_Go_Callback(hObject, eventdata, handles)
% hObject    handle to PB_Go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%0. Create folder to save data in - use selected folder, but make subfolder with current time stamp
save_folder_base = get(handles.ET_SaveFolder,'String');
save_name = get(handles.ET_SaveName,'String');
date_stamp = datestr(datetime,'mm-dd-yyyy_HH-MM-SS');
save_folder = [save_folder_base '\' date_stamp];
handles.save_folder = save_folder;

if not(exist(save_folder,'dir'))
    mkdir(save_folder)
else
    warning('folder already exists')
end

%Optionally start a log file
global DoLog d_num
d_num = 1;
%Possibly start new log file
if DoLog
    logname=[save_folder '\Log_Bayes_' num2str(d_num) '__' datestr(now,'mm_dd_yyyy__HH_MM') '.txt'];
    diary(logname);
end

%Optionally show troubleshooting message on screen
global DoTroubleshoot
if DoTroubleshoot
   disp('Run: PB_Go_Callback...')
end  

global DoExisting %We are doing on demand (use of go button)
DoExisting=0;

%RUNS WHEN GO IS PRESSED
%1. Reset data aquisition
daqreset
if license('test','image_acquisition_toolbox')
    imaqreset
end

rng(17);

%% 

global q
p = gcp();

%% 
restart_mat_paths = get(handles.ET_path_for_restart,'String'); % this could be multiple mats
restart_mat_paths_cell = strsplit(restart_mat_paths,';');

global prior_Seizure_Count prior_Seizure_Duration prior_stim_freq prior_stim_amp prior_stim_width prior_stim_ampRatio prior_stim_chgRatio

%2. Optionally load in data from given location prior to starting current session 
%   Loading in one or more summary files and aggregating data from those.
if not(isempty(restart_mat_paths))

    %Initialize variables
    n_ch_out = length(str2num(get(handles.ET_n_ch_out,'String')));
    prior_Seizure_Count = zeros(1,n_ch_out);
    prior_Seizure_Duration = [];
    prior_stim_freq = [];
    prior_stim_amp = [];
    prior_stim_width = [];
    prior_stim_ampRatio = [];
    prior_stim_chgRatio = [];

    %Loop through previous summary files
    for i_mat = 1:length(restart_mat_paths_cell)
        
        %Access summary file
        mf_restart = matfile(restart_mat_paths_cell{1,i_mat});

        %???!!!???might work with my new saving method
        % prior_Seizure_Count = mf_restart.Seizure_Count; % doesn't work reliably because it is only saved at the end
%         Seizure_Count = min([sum(mf_restart.stim_freq~=0,1); sum(mf_restart.Seizure_Duration~=0,1)],[],1);
        Seizure_Count_mf = mf_restart.Seizure_Count;

        %Loop through channels
        for i_ch = 1:n_ch_out
            if Seizure_Count_mf(i_ch)>0
                %Extract stimulation parameters for current summary file
                Seizure_Duration_ch = mf_restart.Seizure_Duration(1:Seizure_Count_mf(i_ch),i_ch);
                stim_freq_ch = mf_restart.stim_freq(1:Seizure_Count_mf(i_ch),i_ch);
                stim_amp_ch = mf_restart.stim_amp(1:Seizure_Count_mf(i_ch),i_ch);   %tom:was not correctly saved previously due to typo, but data is also in duration amp freq
    %             stim_amp_ch = mf_restart.duration_amp_freq(1:Seizure_Count_mf(i_ch),(i_ch-1)*3+2);
                stim_width_ch = mf_restart.stim_width(1:Seizure_Count_mf(i_ch),i_ch);
                stim_ampRatio_ch = mf_restart.stim_ampRatio(1:Seizure_Count_mf(i_ch),i_ch);
                stim_cghRatio_ch = mf_restart.stim_chgRatio(1:Seizure_Count_mf(i_ch),i_ch);

                %Aggregate stimulation parameters for ALL summary files
                prior_Seizure_Duration(prior_Seizure_Count(i_ch)+(1:Seizure_Count_mf(i_ch)),i_ch) = Seizure_Duration_ch;
                prior_stim_freq(prior_Seizure_Count(i_ch)+(1:Seizure_Count_mf(i_ch)),i_ch) = stim_freq_ch;
                prior_stim_amp(prior_Seizure_Count(i_ch)+(1:Seizure_Count_mf(i_ch)),i_ch) = stim_amp_ch;
                prior_stim_width(prior_Seizure_Count(i_ch)+(1:Seizure_Count_mf(i_ch)),i_ch) = stim_width_ch;
                prior_stim_ampRatio(prior_Seizure_Count(i_ch)+(1:Seizure_Count_mf(i_ch)),i_ch) = stim_ampRatio_ch;
                prior_stim_chgRatio(prior_Seizure_Count(i_ch)+(1:Seizure_Count_mf(i_ch)),i_ch) = stim_cghRatio_ch;
            end
        end
        %Update counts for loading in next file (next iteration)
        prior_Seizure_Count = prior_Seizure_Count+Seizure_Count_mf;
    end

    disp('loaded prior data!')
    
    ET_SF = get(handles.ET_SaveFolder,'String');
    set(handles.ET_SaveFolder,'String',[ET_SF '_RS']);

else
    disp('no prior data loaded')
    prior_Seizure_Count = zeros(1,4);
    prior_Seizure_Duration = [];
    prior_stim_freq = [];
    prior_stim_amp = [];
    prior_stim_width = [];
    prior_stim_ampRatio = [];
    prior_stim_chgRatio = [];
end

%% load hierarchical model
%3. Loading hierarchical model if selected
optim_string = get(handles.ET_bayes_or_AB,'String');

global hier_model_settings P_t_g theta_t_g off_set_est

if strcmpi(optim_string,'Hier') == 1
    
    hier_mat_path = get(handles.ET_hier_model_mat,'String');
    load(hier_mat_path,'theta_t','P_t','n_setting_combos','n_prior_mice','bins_of_optimized_dims','sigma')
    theta_0 = theta_t;
    P_0 = P_t;
    hier_model_settings.n_setting_combos = n_setting_combos;
    hier_model_settings.n_prior_mice = n_prior_mice;
    hier_model_settings.bins_of_optimized_dims = bins_of_optimized_dims;
    hier_model_settings.sigma = sigma;
    clear theta_t P_t
    
    P_t_g = repmat(P_0,1,1,n_ch_out);       % replicate P_0 for every mouse
    theta_t_g = repmat(theta_0,1,n_ch_out); % theta_t is column vector, replicate it into a matrix for each mouse

    global P_orig theta_orig                % make a copy of the initial P and theta so that recalculation can be done from scratch if offset changes
    P_orig = P_t_g;
    theta_orig = theta_t_g;

    off_set_est = log(2)*ones(1,n_ch_out);
end

n_ch_out = length(str2num(get(handles.ET_n_ch_out,'String')));


%% create save folder

save_settings(save_folder, 'settings',handles);

n_cams = str2num(get(handles.ET_n_cams,'String'));
if n_cams >0 %launch new matlab to record videos
    eval(['!matlab /r n_cams=' num2str(n_cams) ';path=''' save_folder '\' save_name '_vid'';cam_corder&'])
end

% eval(['!matlab /r bo_res_path_base=''' save_folder ''';&']);% ND_load_n_blot_bo_SEP_res& % commenting out for debug

%% setup daq
% 'PCI-6251'
% 'USB-6229 (BNC)'
% 'USB-6221 (BNC)' ?
devices = daq.getDevices;
device_ID_str = get(handles.ET_device_ID,'String');
device_ID_stim = get(handles.ET_stim_dev,'String');

for i = 1:length(devices)
    if strcmp(devices(i).Model, device_ID_str)
        NI_dev = devices(i);
    end
    if strcmp(devices(i).Model,device_ID_stim)
        stim_dev = devices(i);
    end
end

s = daq.createSession('ni');
% s_stim = daq.createSession('ni');

handles.s = s;
% handles.s_stim = s_stim;

global fs_rec 
fs_rec = str2num(get(handles.ET_fs,'String'));
% fs_stim = str2num(get(handles.ET_fs_stim,'String'));

s.Rate = fs_rec;
s.IsContinuous = true;

% s_stim.Rate = fs_stim;
% s_stim.IsContinuous = true;

ch_out_vec = str2num(get(handles.ET_n_ch_out,'String'));
n_ch_out = length(ch_out_vec);

for i_ch = 1:n_ch_out
    q{1,i_ch} = parallel.pool.PollableDataQueue;
end

% % stim_on_durs = randi(35,10,1);
% % stim_off_durs = randi(35,10,1);
% % for i_ch_out = 1:n_ch_out
% %     figure(i_ch_out)
% %     subplot(1,2,1)
% %     histogram(stim_on_durs,.5:1:100)
% %     title('stim')
% %     subplot(1,2,2)
% %     histogram(stim_off_durs,.5:1:100)
% %     title('no stim')
% % end

global stim_remaining
global stim_flag

seizure_detection_ch_vec = str2num(get(handles.ET_seizure_detection_channels,'String'));% detect seizures on ACH0 and ACH2, these are BNC-2090 numbers
  
if length(seizure_detection_ch_vec) ~= n_ch_out
    error('length(seizure_detection_ch_vec) ~= n_ch_out')
end

addAnalogOutputChannel(s, stim_dev.ID, ch_out_vec, 'Voltage');
  
for i_D_out = 1:length(ch_out_vec)
    addDigitalChannel(s,stim_dev.ID,['port0/line' num2str(i_D_out-1)],'OutputOnly');
end

ch_in_vec = str2num(get(handles.ET_n_ch_in,'String')); % hard coded 1 input channel
n_ch_in = length(ch_in_vec);

for i_ch = ch_in_vec % add ephys recording channels
    addAnalogInputChannel(s,NI_dev.ID, i_ch, 'Voltage');
end

if or(strcmpi(device_ID_str,'USB-6343 (BNC)'),strcmpi(device_ID_str,'USB-6229 (BNC)'))
    for i_D_out = 1:length(ch_out_vec)
        addDigitalChannel(s,NI_dev.ID,['port0/line' num2str(i_D_out-1+10)],'InputOnly');
    end
elseif or(strcmpi(device_ID_str,'USB-6221 (BNC)'),strcmpi(device_ID_str,'USB-6341 (BNC)'))
    for i_D_out = 1:length(ch_out_vec)
        addDigitalChannel(s,NI_dev.ID,['port0/line' num2str(i_D_out-1+4)],'InputOnly'); %6341
    end
end

%%
%Initialize data structures needed
InitializeVarsInput(n_ch_in,n_ch_out);        %Initializes variables used in spike detection and seizure detection

%call global variables that will be needed in the rest of this function
global fast_int slow_int Seizure_On Seizure_Off Seizure_Duration
global spike_count
global spikes_trigger spikes_trigger_stop
global last_spike_time Seizure_Start_First_Spike_Time
global time_out
global in_chunk out_chunk
global ro_data co_data ro_dd co_dd ro_ard co_ard    %keep track of row number in output data var
global datatype
global TimeStamp_postfix

fast_slow_ratio_trigger = zeros(1,n_ch_out,'logical');

global next_freq next_amp next_width  next_ampRatio next_chgRatio
next_freq = 10*ones(1,n_ch_out); % initialize arbitrarily to 10 Hz
next_amp = 5*ones(1,n_ch_out); % initialize arbitrarily to 1 unit
next_width = 0.150*ones(1,n_ch_out); % initialize arbitrarily to 1 unit
next_ampRatio = 1*ones(1,n_ch_out); % initialize to a 4:1 amplitude ratio
next_chgRatio = 1*ones(1,n_ch_out); % initialize to a 1:1 chg ratio

set(handles.T_NextFreq,'String',num2str(next_freq));
set(handles.T_NextAmp,'String',num2str(next_amp));
set(handles.T_NextWidth,'String',num2str(next_width));

%% initialize matfile
%Data variables get initialized to zeros and saved to file
time=zeros(1,1);                                 %EFF21: time as separate variable
data=zeros(1,n_ch_in+n_ch_out);                  %EFF21: next n_ch_in columns are the input channels, and the last n_ch_out columns are the output channels
art_rem_data = zeros(1,n_ch_out,datatype);       %EFF21
art_rem_data_raw = zeros(1,n_ch_out,datatype);   %EFF21
detect_data = zeros(1,n_ch_out,datatype);    %EFF21
detect_data_raw = zeros(1,n_ch_out,datatype);%EFF21
ro_data=1; co_data=size(data,2);                        %Keep track of row we are in
ro_dd=1; co_dd=size(detect_data,2);                     %Keep track of row we are in
ro_ard=1; co_ard=size(art_rem_data,2);             %Keep track of row we are in

%Define extra variables to be stored in the files, to be used for
%identification purposes.
LiveRecording=1;                                          %LiveRecording=0 means offline reanalysis, LiveRecording=1 means new recording
ReferenceDir='';                                          %This means the _d_ analysis files are in this same directory
ID=save_name;                                             %identifying ID (???for now this is the save_name chosen
SettingsName='settings.mat';                              %name of the settings file to be used
TimeStampString=datestr(datetime,'mm-dd-yyyy_HH-MM-SS');  %current time stamps. Included in EVERY file saved.
TimeStamp=now;                                            %Current time and date
TimeStamp_postfix='';                                     %[' (' TimeStamp ')'];

%Define extra variables so that individual files can store some of the data
%as a safety precaution when the large file gets corrupted.
global last_spike_time_Current Seizure_Start_First_Spike_Time_Current Seizure_Duration_Current stim_Current
last_spike_time_Current =  zeros(1,n_ch_out);
Seizure_Start_First_Spike_Time_Current = zeros(1,n_ch_out);
Seizure_Duration_Current = zeros(1,n_ch_out);
stim_Current = zeros(1,5*n_ch_out); %n_ch_out*[freq,amp,width,ampratio,chgratio]

%Create the output files and preload them with variabes to be adjusted
%throughout.
global mf mf2 save_mat_path2_base save_mat_path save_mat_path2
save_mat_path = [save_folder '\' save_name '.mat']
save_mat_path2_base = [save_folder '\' save_name '_d_'] % c
save_mat_path2 = [save_mat_path2_base  num2str(d_num) TimeStamp_postfix '.mat'];

save(save_mat_path2,'time','data','art_rem_data','art_rem_data_raw','detect_data','detect_data_raw','fs_rec',...
                    'ID','TimeStamp','last_spike_time_Current','Seizure_Start_First_Spike_Time_Current','Seizure_Duration_Current','stim_Current','-v7.3','-nocompression')
save(save_mat_path,'fast_int','slow_int','Seizure_Off','stim_flag','Seizure_Duration','spike_count','fast_slow_ratio_trigger','spikes_trigger','spikes_trigger_stop','Seizure_Start_First_Spike_Time','last_spike_time','time_out',...
                   'LiveRecording','ReferenceDir','ID','SettingsName','TimeStamp','TimeStampString','-v7.3','-nocompression')

%EFF21: No need anymore since I rewrote structure
%spike_count_history = spike_count_history(:,:,1); % remove 3rd dimension
%spike_count_history_raw = spike_count_history_raw(:,:,1); % remove 3rd dimension
%spike_count_history_combined = spike_count_history_combined(:,:,1); % remove 3rd dimension
%We initialized these vars to large size for saving to file, but now
%these can be reduced for internal use, where we only need to save most
%recent 10 rows. 11 is HARDCODED as 1+10
clear data spike_count

mf = matfile(save_mat_path,'Writable',true);
mf2 = matfile(save_mat_path2,'Writable',true);

%% specify optimized variables
% opt_freq = []; % Hz
% opt_amp = []; % Volts?


%% (don't) run the bayes opt and plots one time to get everything compiled, just set the first point


% freqLo_vec = str2num(get(handles.ET_freqLo,'String'));
% freqHi_vec = str2num(get(handles.ET_freqHi,'String'));
% ampLo_vec = str2num(get(handles.ET_ampLo,'String'));
% ampHi_vec = str2num(get(handles.ET_ampHi,'String'));
% % widthLo_vec = str2num(get(handles.ET_widthLo,'String'));
% % widthHi_vec = str2num(get(handles.ET_widthHi,'String'));
% width_vec = str2num(get(handles.ET_widthLo,'String'));
% 
% InitialObjective = [15 31 5 52]'; % make up some data
% InitialObjective = log(InitialObjective);
% freq = randi([freqLo_vec(1), freqHi_vec(1)],4,1);
% amp = randi(fix([ampLo_vec(1), ampHi_vec(1)]),4,1);;
% InitialX = table(freq, amp);
% 
% for i_ch = 1:n_ch_out
%     amp_range = [ampLo_vec(i_ch) ampHi_vec(i_ch)];
%     freq_range = [freqLo_vec(i_ch) freqHi_vec(i_ch)];
%     opt_amp = Variable('amplitude',amp_range,'Type','integer'); % Volts
%     opt_freq = optimizableVariable('frequency',freq_range); % Hz
%     parfeval(@BO_wrapper,0,opt_freq, opt_amp, InitialX, InitialObjective, q{1,i_ch});
% end
% 
% 
% for i_ch_out = 1:n_ch_out
%     gotMsg = 0;
%     while gotMsg ~= 1
%         pause(.1)
%         [res, gotMsg] = poll(q{1,i_ch_out}, .02); % should save each res
% 
%         if gotMsg
%     %         close all
% % %             tic  ;% plots are now handled seperately
% % %             figure(i_ch_out)
% % %     %         plot(res,@plotObjectiveModel) %  A_BPO_vis, would be good to make these into a video
% % %             amp_range = [ampLo_vec(i_ch) ampHi_vec(i_ch)];
% % %             freq_range = [freqLo_vec(i_ch) freqHi_vec(i_ch)];
% % %             plot_bo(res, freq_range, amp_range, [0 4.5], 35)
% % %             toc
% 
%     %         eval(['res_ch_' num2str(i_ch_out) '_sz_' num2str(Seizure_Count(1,i_ch_out)) '= res;'])
%     %         tic
%     %         save(save_mat_path,['res_ch_' num2str(i_ch_out) '_sz_' num2str(Seizure_Count(1,i_ch_out))],'-nocompression','-append')
%     %         toc
% 
%             next_freq(1,i_ch_out) = res.NextPoint{1,1};
%             next_amp(1,i_ch_out) = res.NextPoint{1,2};
%             next_width(1,i_ch_out) = res.NextPoint{1,3};
% 
%             set(handles.T_NextFreq,'String',num2str(next_freq));
%             set(handles.T_NextAmp,'String',num2str(next_amp));
%             set(handles.T_NextWidth,'String',num2str(next_width));
%         end
%     end
% end

next_freq(1,:) = 10;
next_amp(1,:) = 5;
next_width(1,:) = 0.1;
next_ampRatio(1,:) = 1;
next_chgRatio(1,:) = 1;

set(handles.T_NextFreq,'String',num2str(next_freq));
set(handles.T_NextAmp,'String',num2str(next_amp));
set(handles.T_NextWidth,'String',num2str(next_width));
set(handles.T_NextAmpRatio,'String',num2str(next_ampRatio));
set(handles.T_NextChgRatio,'String',num2str(next_chgRatio));

%% add listeners
lh1 = addlistener(s,'DataAvailable', @(src, event) Process_Plot_Save(src,event, handles.A_MainPlot,ch_in_vec,seizure_detection_ch_vec,n_ch_in,n_ch_out,save_mat_path, handles) );
lh2 = addlistener(s,'DataRequired', @(src, event) Generate_Stim_Vec(src,event,handles));
lh3 = addlistener(s, 'ErrorOccurred', @(~,event) disp(getReport(event.Error)));

s.NotifyWhenDataAvailableExceeds = fix(fs_rec*in_chunk); % input buffer treshold, hard coded
%CKM: When there is 1s (in_chunk) of data available, use listener lh1 to run Process_Plot

s.NotifyWhenScansQueuedBelow = fix(fs_rec*1.4*out_chunk); % output buffer threshold, hard coded
%CKM: When there is less than 1.4s of data available in the buffer, it will trigger Generate_Stim_Vec using listener lh2

%% buffer 5 seconds of zeros to the daq output to start with, then wait 0.5second
stim_vec_init = zeros(fix(5*fs_rec),n_ch_out*2);
queueOutputData(s, [stim_vec_init]);
pause(0.5);

%% start daq
disp('starting daq')
daq_start_now = now;
daq_start_datetime = datetime;

% startBackground(s_stim);
startBackground(s); % start the daq

%% record daq start time info
[Y,M,D,H,MN,S] = datevec(daq_start_datetime);
daq_date_vector = [Y,M,D,H,MN,S];
save(save_mat_path,'daq_date_vector','daq_start_now','daq_start_datetime','-nocompression','-append')

% handles.vid = vid;
% handles.n_cams = n_cams;

guidata(hObject,handles);


%% --- Executes on button press in PB_Stop: PRESSING THE STOP BUTTON
function PB_Stop_Callback(hObject, eventdata, handles)
% hObject    handle to PB_Stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DoSaveToFile

global DoTroubleshoot
if DoTroubleshoot
   disp('PB_Stop_Callback...')
end

%CKM: Called when STOP button is pressed.
s = handles.s;
s.stop()
% s_stim = handles.s_stim;
% s_stim.stop()

ch_out_vec = str2num(get(handles.ET_n_ch_out,'String'));
n_ch_out = length(ch_out_vec);

if DoSaveToFile
    %Final save to file of all vars that did not yet get updated since last
    %switch of data-files.
    UpdateFiles(n_ch_out,1);
end

% for i_ch_out = 1:n_ch_out
%     saveas(i_ch_out,[handles.save_folder '\figure_' num2str(i_ch_out) '.fig'])
%     saveas(i_ch_out,[handles.save_folder '\figure_' num2str(i_ch_out) '.png'])
% end

%Turn off diary if it was running
global DoLog
if DoLog
    diary off
end


%% --- Executes when you want to open existing file to check for seizures under current settings
function btn_open_Callback(hObject, eventdata, handles)
%When using this option:
%- data is loaded from file, not from daq
%- it does not save to file any of your data (might later)
global DoStopOffline DoRestartOffline
global DoExisting   %1=Load from existing (offline), 0=online
global DoSaveToFile %0=don't save variables to file, 1=save to file
global fs_rec       %sampling rate
DoExisting = 1;     %1=we do offline run on existing file
DoStopOffline = 0;  %1=stop offline analysis after current chunk
DoRestartOffline = 0;%0=load in data for first time, 1=restart same data
DoSaveToFile = 0;

%Ask user to select one or more files
global ReanalysisFolder ReanalysisFiles ReanalysisIndex
[ReanalysisFiles,ReanalysisFolder]=uigetfile('*.mat','MultiSelect','on');
if ~iscell(ReanalysisFiles)
    if ReanalysisFiles==0
        return
    else
        ReanalysisFiles={ReanalysisFiles};
        ReanalysisIndex=1;
    end
end
ReanalysisPath = [ReanalysisFolder ReanalysisFiles{1}];

%Extract necessary input from file1 to obtain some basic info (such as fs)
InputMatfile=matfile(ReanalysisPath);
%Check if time exists as separate variable or is included as first column of data
listOfVariables = who('-file', ReanalysisPath);
TimeIncluded=0+~ismember('time', listOfVariables); % returns true if it is first column

SizeData=size(InputMatfile,'data');
SizeNrChOut=size(InputMatfile,'detect_data');
SizeNrChOut=SizeNrChOut(2)-TimeIncluded;            %nr output channels
SizeNrChIn=SizeData(2)-SizeNrChOut-TimeIncluded;    %nr input channels
fs_rec=InputMatfile.fs_rec;                         %sampling rate

%Unhide popup
set(handles.PanelTest,'Visible','on');

%Update parts of GUI to reflect we loaded in data with certain criteria
set(handles.ListOfFiles,'String',ReanalysisFiles);
set(handles.ListOfFiles,'value',1);
set(handles.TextFS,'String',['FS: ' num2str(fs_rec) ' Hz']);
locs=strfind(ReanalysisFolder,'\');
set(handles.TextFN,'String',['Loaded: ' ReanalysisFolder(1:3) '...' ReanalysisFolder(locs(end-1):end) ReanalysisFiles{ReanalysisIndex}]);
set(handles.TextChIn,'String',['#Channels In: ' num2str(SizeNrChIn)]);
set(handles.TextChOut,'String',['#Channels Out: ' num2str(SizeNrChOut)]);
set(handles.InputFileName,'String',[ReanalysisFolder 'Reanalysis']);

%Check if the number of channels specified in settings matches what we
%pulled from the files. if not, show warning
numberin = str2num(get(handles.ET_n_ch_in,'String'));% detect seizures on ACH0 and ACH2, these are BNC-2090 numbers
if length(numberin)~=SizeNrChIn
    warndlg('The number of input channels in datafile does not match "ch in vec" in the GUI')
end

n_ch_out = length(str2num(get(handles.ET_n_ch_out,'String')));
if n_ch_out~=SizeNrChOut
    warndlg('The number of output channels in datafile does not match "ch out vec" in the GUI')
end

%Disable open button
set(handles.btn_open,'Enable','off');
set(handles.ButtonStart,'Enable','on');

%reset figure
hold(handles.A_MainPlot, 'off')
plot(handles.A_MainPlot,0,0);
hold(handles.A_MainPlot, 'on')


%% --- Initialize variables (called by both online and offline detections)
function InitializeVarsInput(n_ch_in,n_ch_out)
% Contains variables needed for detection of spikes and events
global fs_rec   %recording frequency (sampling rate)
global n_read   %number of iterations (chunks)
n_read=1;

% Setup how many times the sampling rate it loaded in at a time and set up at a time
global out_chunk in_chunk
out_chunk = 1; % these are very long, 1/2 second would be better, but computation time of bayes opt is too slow to fit inside the loop
in_chunk = 1;

global fast_int slow_int 
global Seizure_On Seizure_Off Seizure_Count Seizure_Duration Seizure_Start_Ind 
global Seizure_Count_Current
global duration_amp_freq
global stim_flag stim_remaining
global Sz_On_Index Sz_Off_Index Sz_Index
global datatype

datatype='single';                                   %EFF21: save all data files except for the raw data as single format
fast_int = zeros(1,n_ch_out);
slow_int = zeros(1,n_ch_out);
Seizure_On = false(2,n_ch_out);                      %EFF21: Only keep track of last 2 entries, don't save to file
Seizure_Off = false(1,n_ch_out);                     %EFF21: In file, will be sparse logical, only updated once every file
stim_flag = false(1, n_ch_out);                      %EFF21: In file, will be sparse logical, only updated once every file
Sz_On_Index = [];                                    %EFF21: (row,col) of ON, currently not prefilled to length (used for stim_flag creation)
Sz_Off_Index = [];                                   %EFF21: (row,col) of OFF, currently not prefilled to length (used for Seizure_Off creation)
Sz_Index = [1,1];                                    %EFF21: index where current file starts for on and off
Seizure_Count = zeros(1,n_ch_out,'single');          %keep track of #seizures in ALL files %EFF21
Seizure_Count_Current = zeros(1,n_ch_out,'single');  %keep track of #seizures in current file %EFF21
Seizure_Duration = zeros(1,n_ch_out);
Seizure_Start_Ind = zeros(1,n_ch_out);
duration_amp_freq = zeros(1,3*n_ch_out);
stim_remaining = zeros(1,n_ch_out);

maxsize=60*10*200;    %400 files initialization

%Set up spike detection parameters
global pos_spike_count neg_spike_count pos_spike_count_raw neg_spike_count_raw valid_spike_count
global spike_count_history spike_count_history_raw spike_count_history_combined spike_count
pos_spike_count = zeros(1,n_ch_out,'uint32');
neg_spike_count = zeros(1,n_ch_out,'uint32');
pos_spike_count_raw = zeros(1,n_ch_out,'uint32');
neg_spike_count_raw = zeros(1,n_ch_out,'uint32');
valid_spike_count = zeros(1,n_ch_out,'uint32');
%HARDCODED: if in_chunck ever changes to value ~=1s, then following will need to be changed.
%Also, indexing tied to initialize to 10chunks, is hardcoded elsewhere.
spike_count = zeros(maxsize,3*n_ch_out,'uint16');            % save spike counts into one large array in file to reduce #saves to file
spike_count_history = zeros(10,n_ch_out,'uint16');           % spike_count_history is initialized to 10 seconds, extra dimension to make file 3D
spike_count_history_raw = zeros(10,n_ch_out,'uint16');       % spike_count_history is initialized to 10 seconds, extra dimension to make file 3D
spike_count_history_combined = zeros(10,n_ch_out,'uint16');  % spike_count_history is initialized to 10 seconds, extra dimension to make file 3D
%EFF21:Changing to 2D array from 3D array
global spikes_trigger spikes_trigger_stop spike_times_buffer_cell 
global Seizure_Start_First_Spike_Time last_spike_time time_out
%TODO EFF21: use sparse creation periodically
spikes_trigger = false(10,n_ch_out);        %EFF21 change to 10dim array. Will use circshift and always change 10th entry
spikes_trigger_stop = false(10,n_ch_out);   %EFF21 change to 10dim array. Will use circshift and always change 10th entry
global Tr_On_Index Tr_Off_Index Tr_Index
Tr_On_Index =[];            %EFF21: 
Tr_Off_Index=[];            %EFF21: 
Tr_Index=[1,1];             %EFF21: 
spike_times_buffer_cell = cell(5,n_ch_out); % 2 seconds hard coded, most recent spike times from teh past two seconds, newest goes into end, n_ch_out
last_spike_time =  zeros(1,n_ch_out);
Seizure_Start_First_Spike_Time = zeros(1,n_ch_out);
time_out = 10*ones(1,n_ch_out,'int32')+int32((1:n_ch_out)*10); % added some stagger for the start %EFF21 NO NEED TO SAVE TO FILE

global ready_for_next_sz
global prev_digi_data prev_digi_time
global previous_data_chunk
global prev_data_chunk_len 
ready_for_next_sz = true(1,n_ch_out);
prev_digi_data = zeros(100,n_ch_out);
prev_digi_time = (-100:-1)'/fs_rec;
prev_data_chunk_len = fix(fs_rec*0.05);
previous_data_chunk = zeros(prev_data_chunk_len,n_ch_in+n_ch_out);

%Set up filters
global adaptive_filt_cell adaptive_filt_cell_raw
adaptive_filt_cell=cell(n_ch_out);            %preallocate cell array
adaptive_filt_cell{n_ch_out,n_ch_out}=[];     %preallocate cell array
adaptive_filt_cell_raw=cell(n_ch_out);        %preallocate cell array
adaptive_filt_cell_raw{n_ch_out,n_ch_out}=[]; %preallocate cell array
LMS_filter_len = fix(fs_rec*0.004);
LMS_filter_len_raw = fix(fs_rec*0.025);
for i_ch = 1:n_ch_out
    for j_ch = 1:n_ch_out
        adaptive_filt_cell{i_ch,j_ch}.filter = dsp.LMSFilter('Length',LMS_filter_len,'StepSizeSource','Input port','AdaptInputPort',true,'WeightsResetInputPort',true);
        adaptive_filt_cell_raw{i_ch,j_ch}.filter = dsp.LMSFilter('Length',LMS_filter_len_raw,'StepSizeSource','Input port','AdaptInputPort',true,'WeightsResetInputPort',true);
    end
end


%% --- Executes using listener lh1 whenever 1s of new data is available in the buffer or in offline detection from file
function Process_Plot_Save(src,event,plot_handle, ch_in_vec, seizure_detection_ch_vec, n_ch_in, n_ch_out, save_mat_path, handles) 
%This function can be called during online detection or offline detection.
%When it is called in offline mode, it will not save anything to file and it will obtain all data from loaded file
%It will also not create any output pulses in this mode

%Optionally show troubleshooting code in command window
global DoTroubleshoot
if DoTroubleshoot
   disp('Run: Process_Plot_Save (Called by listener lh1)...')
end
disp('..........')

%Check if this is online detection or offline detection
global DoExisting %1=offline, 0=online
global DoSaveToFile
if ~DoExisting
    DoSaveToFile=1; %ALWAYS SAVE TO FILE WHEN REAL TIME
    set(handles.CheckSave,'Value',1);
end

global DoResaveDataFiles DoSaveData DoFinishedFile
%if DoExisting=1, but we want to resave data files ANYWAY, then we use
%"DoResaveDataFiles". So then ~DoExisting || DoResaveDataFiles
DoSaveData=(~DoExisting || DoResaveDataFiles);

%Name all global variables to be loaded in for this function
global stim_flag
global fs_rec
global duration_amp_freq
global q
global n_read fast_int slow_int Seizure_On Seizure_Off Seizure_Count Seizure_Duration Seizure_Start_Ind 
global Seizure_Count_Current
global next_freq next_amp next_width next_ampRatio next_chgRatio
global in_chunk spike_count_history spike_count_history_raw spike_count_history_combined spike_count
global pos_spike_count neg_spike_count pos_spike_count_raw neg_spike_count_raw valid_spike_count
global spike_times_buffer_cell
global Seizure_Start_First_Spike_Time last_spike_time
global time_out
global spikes_trigger spikes_trigger_stop
global prev_filt_data_x prev_filt_data_y
global adaptive_filt_cell adaptive_filt_cell_raw
global stim_amp stim_freq stim_width stim_ampRatio stim_chgRatio
global ready_for_next_sz
global last_spike_time_Current Seizure_Start_First_Spike_Time_Current Seizure_Duration_Current stim_Current
global ro_data co_data ro_dd co_dd ro_ard co_ard
global Sz_On_Index Sz_Off_Index Sz_Index
global Tr_On_Index Tr_Off_Index Tr_Index
global datatype
global prior_Seizure_Count prior_Seizure_Duration prior_stim_freq prior_stim_amp prior_stim_width prior_stim_ampRatio prior_stim_chgRatio
global save_mat_path2
global TimeStamp_postfix
n_read = n_read+1; %index of calls to Process_Plot_Save. Gets initialized at 1 when pressing GO, so has value 2 first time this is run

%% start new matlab data file 'mf2' every so often (Data_d_1, Data_d_2,...)
global d_num mf mf2 save_mat_path2_base DoLog
ch_out_vec = str2num(get(handles.ET_n_ch_out,'String'));
n_ch_out = length(ch_out_vec);

%If in online mode, setup the next file here (this fn is called via listener. In offline mode we automatically
%repeat until whole file is read in function ButtonStart and call this
%function in a loop from there.
%1. FINISH UP PREVIOUS FILE - bookkeeping: call UpdateFiles to save final
%vars for previous read to files (both summary and data files).
%REAL TIME: mod(n)read,60*10 && ~DoExisting => do this every 600 chunks of data (every 600s)
%OFFLINE: => do this when a file is completely read and finished, 
if (mod(n_read,60*10)==0 && ~DoExisting) || (DoFinishedFile && DoExisting)
    %if we are rerunning, we really only want to do this 
    UpdateFiles(n_ch_out,DoSaveData);
end

%2. SET UP NEW FILE
%This happens every 600 runs of Process_Plot_Save which occurs each second of data coming in
if (mod(n_read,60*10)==0 && ~DoExisting) || (DoFinishedFile && (DoExisting && DoResaveDataFiles)) 
    d_num = d_num + 1;                      %Index of the data file. Gets initialized at 1 when pressing GO
    
    %Display message to show where we are at
    %disp(['n_read=' num2str(n_read) ', d_num= ' num2str(d_num)]);
    
    %Possibly start new log file (if DoLog=2)
    if DoLog>1
        diary off;  %stop diary from previous file
        logname=['Log_Bayes_' num2str(d_num) '__' datestr(now,'mm_dd_yyyy__HH_MM') '.txt'];
        diary(logname); %start new diary
    end
    
    %Initialize stat variables for new file
    for i_ch=1:n_ch_out %Reset counts for #seizures for current file. Be careful, because sometimes last spike occurs in next file, so sometimes we should initialize to 1.
        if size(last_spike_time_Current,1)<Seizure_Count_Current(1,i_ch)
            disp(['Carrying over seizure (last spike) started in previous file on channel ' num2str(i_ch)]);
            Seizure_Count_Current(1,i_ch)=1;  %last spike on channel did not yet occur, so start index on 1 for next file 
        elseif Seizure_Count_Current(1,i_ch)==0
            Seizure_Count_Current(1,i_ch)=0; %no seizure activity at all on this channel
        elseif last_spike_time_Current(Seizure_Count_Current(1,i_ch),i_ch)==0
            Seizure_Count_Current(1,i_ch)=1;  %last spike on channel did not yet occur, so start index on 1 for next file 
            disp(['Carrying over seizure (last spike) started in previous file on channel ' num2str(i_ch)]);
        else
            Seizure_Count_Current(1,i_ch)=0;  %last spike already occurred, so reset to 0
        end
    end
    last_spike_time_Current =  zeros(1,n_ch_out);
    Seizure_Start_First_Spike_Time_Current = zeros(1,n_ch_out);
    Seizure_Duration_Current = zeros(1,n_ch_out);
    stim_Current = zeros(1,5*n_ch_out);

    %Initialize data structures for new file
    time=zeros(1,1);
    data=zeros(1,n_ch_in+n_ch_out);            %Column 1 is time stamps, next n_ch_in columns are the input channels, and the last n_ch_out columns are the output channels
    art_rem_data = zeros(1,n_ch_out,datatype);
    art_rem_data_raw = zeros(1,n_ch_out,datatype);
    detect_data = zeros(1,n_ch_out,datatype);
    detect_data_raw = zeros(1,n_ch_out,datatype);
    ro_data=1; co_data=size(data,2);                        %Keep track of row we are in
    ro_dd=1; co_dd=size(detect_data,2);                     %Keep track of row we are in
    ro_ard=1; co_ard=size(art_rem_data,2);             %Keep track of row we are in
    
    ID=get(handles.ET_SaveName,'String');                     %identifying ID (???for now this is the save_name chosen
    TimeStampString=datestr(datetime,'mm-dd-yyyy_HH-MM-SS');  %current time stamps. Included in EVERY file saved.
    TimeStamp=now;
    %TimeStamp_postfix is global, and defined in earlier setup
    save_mat_path2 = [save_mat_path2_base num2str(d_num) TimeStamp_postfix '.mat'] % c
    save(save_mat_path2,'time', 'data','art_rem_data','art_rem_data_raw','detect_data','detect_data_raw','fs_rec',...
                    'ID','TimeStamp','TimeStampString','last_spike_time_Current','Seizure_Start_First_Spike_Time_Current','Seizure_Duration_Current','stim_Current','-v7.3','-nocompression')
    
    clear time data art_rem_data detect_data art_rem_data_raw
    mf2 = matfile(save_mat_path2,'Writable',true);
end


%% save next block of raw data and time
%EFF21 - don't keep requesting dimensions of data, but rather keep track of
%it in variables in workspace.
[r, c] = size(event.Data);
if DoSaveData
    if co_data~=c
        error('data dims not equal in columns')
    end
    mf2.time(ro_data+(1:r), :) = event.TimeStamps;   %EFF21, separate out time
    mf2.data(ro_data+(1:r), :) = event.Data;
    ro_data=ro_data+r;
end

%% remove artifacts
detect_ch_logical = ismember(ch_in_vec,seizure_detection_ch_vec);
detect_ch_index = find(detect_ch_logical);
stim_ch_logical = [0, detect_ch_logical(1:end-1)];
stim_ch_index = find(stim_ch_logical);

art_rem_data = zeros(r,n_ch_out);
art_rem_data_raw = zeros(r,n_ch_out);

n_ch = size(event.Data,2);
digital_ch_index = length(ch_in_vec)+(1:n_ch-length(ch_in_vec));

%TODO EFF21, take butter definition out of loop. Does not have to be redefined EVERY time processplot is called
[b_H, a_H] = butter(1,250/(fs_rec/2),'high');
[b_S, a_S] = butter(2,[55 65]/(fs_rec/2),'stop');
[b_H_raw, a_H_raw] = butter(1,5/(fs_rec/2),'high');
[b_L_raw, a_L_raw] = butter(1,80/(fs_rec/2),'low');

global zf0 zf01 zf0_raw zf01_raw zf02_raw
global previous_data_chunk prev_data_chunk_len

estim_or_opto = get(handles.ET_estim_or_opto,'String');

for i_ch = 1:n_ch_out
if any(isnan(previous_data_chunk))
    error('previous data chunk has NaNs')
end
if any(isnan(event.Data))
    error('incoming data (event.Data) has NaNs')
end
    
    pre_padded_Data = [previous_data_chunk; event.Data];
    
    previous_data_samples = size(previous_data_chunk,1);
    
    data = pre_padded_Data(:,detect_ch_index(i_ch));
    data_raw = data;
    
    
    %% HP filter for HF
    if n_read == 2 %this is first 600s file of data gathered
        [data,zf0] = filter(b_H,a_H,data);
        [data,zf01] = filter(b_S,a_S,data);
    else %future cycles of obtaining data
        [data,zf0] = filter(b_H,a_H,data,zf0);
        [data,zf01] = filter(b_S,a_S,data,zf01);
    end    

    train_a_mat = event.Data(:,digital_ch_index);

    %% adaptive filter
    
    data_bckup = data;
    data_bckup_raw = data_raw;
    
    if or(any(isnan(data_bckup)), any(isnan(data)))
        error('data or backup contains nans')
    end
    if or(any(isnan(data_bckup_raw)), any(isnan(data_raw)))
        error('data_raw or raw backup contains nans')
    end
    
    for j_ch = 1:n_ch_out  % should this be in a for loop????
        %% load reference signal and train signal
        ref = pre_padded_Data(:,stim_ch_index(j_ch)); % prepend a bit of the previous second's data
        ref(abs(ref)<0.02) = 0; % maybe should measure and remove the DC instead to be more linear

        a_logical_vec = mean(train_a_mat,1)>0.2;

        a_logical = and(a_logical_vec(j_ch), not(any(a_logical_vec(setdiff(1:n_ch_out,j_ch)))));

        %% LMS
        
        if strcmpi(estim_or_opto,'estim')
            StepSize = 0.005; %previously was 0.01, decreased to promote stability
        elseif strcmpi(estim_or_opto,'opto')
            StepSize = 0.0000001; % very small step size if opto to promote stability, it would be better to disable artifact removal for opto
        else
            error('must be estim or opto')
        end

        reset_flag = false;
        [~, ~, ~] = adaptive_filt_cell{i_ch,j_ch}.filter(ref,data,StepSize,a_logical,reset_flag); % run it once and toss it to improve model  
        [~, ~, ~] = adaptive_filt_cell{i_ch,j_ch}.filter(ref,data,StepSize,a_logical,reset_flag); % run it twice and toss it to improve model   
        [~, data, ~] = adaptive_filt_cell{i_ch,j_ch}.filter(ref,data,StepSize,a_logical,reset_flag);
        
        StepSize_raw = 0.00175; %decreased from 0.0035 to promote stability
        [~, ~, ~] = adaptive_filt_cell_raw{i_ch,j_ch}.filter(ref,data_raw,StepSize_raw,a_logical,reset_flag); % run it once and toss it to improve model   
        [~, ~, ~] = adaptive_filt_cell_raw{i_ch,j_ch}.filter(ref,data_raw,StepSize_raw,a_logical,reset_flag); % run it twice and toss it to improve model   
        [~, data_raw, ~] = adaptive_filt_cell_raw{i_ch,j_ch}.filter(ref,data_raw,StepSize_raw,a_logical,reset_flag);

        nan_loop_ctr = 0;
        while any(isnan(data))
            disp('NaNs!! resetting')
            
            if any(isnan(data_bckup))
                error('backup contains nans too, so this can''t be fixed')
                %  could maybe set all nans to zero and proceed
            end
            
            if any(isnan(ref))
                error('ref contains nans, so this can''t be fixed')
                %  could maybe set all nans to zero and proceed
            end
            
            nan_loop_ctr = nan_loop_ctr + 1
            reset_flag = true;
            
%             if any(isnan(data(1:previous_data_samples,:)))
%                 warning('nan in prev data, zeroing')
%                 ref(1:previous_data_samples,:) = 0; % must remove the previous padding bc it contains nans
%                 data(1:previous_data_samples,:) = 0;
%             end
            
            [~, ~, ~] = adaptive_filt_cell{i_ch,j_ch}.filter(ref,data_bckup,0.05^nan_loop_ctr*StepSize,a_logical,reset_flag); % start over from backup copy
            reset_flag = false;
            [~, ~, ~] = adaptive_filt_cell{i_ch,j_ch}.filter(ref,data_bckup,0.05^nan_loop_ctr*StepSize,a_logical,reset_flag); % start over from backup copy
            [~, data, ~] = adaptive_filt_cell{i_ch,j_ch}.filter(ref,data_bckup,0.05^nan_loop_ctr*StepSize,a_logical,reset_flag); % start over from backup copy
            
            if nan_loop_ctr > 10
                error('cant fix nans')
            end
       
        end    
        
        data_bckup = data; % it made it through, so update the backup
        
        nan_loop_ctr = 0;
        while any(isnan(data_raw))
            disp('NaNs in raw!! resetting')
            
            if any(isnan(data_bckup_raw))
                error('backup contains nans too, so this can''t be fixed')
            end
            
            if any(isnan(ref))
                error('ref contains nans, so this can''t be fixed')
                %  could maybe set all nans to zero and proceed
            end
             
            nan_loop_ctr = nan_loop_ctr + 1
            reset_flag = true;
            
%             if any(isnan(data_raw(1:previous_data_samples,:)))
%                 warning('nan in prev data, zeroing')
%                 ref(1:previous_data_samples,:) = 0; % must remove the previous padding bc it contains nans
%                 data_raw(1:previous_data_samples,:) = 0;
%             end
            
            [~, ~, ~] = adaptive_filt_cell_raw{i_ch,j_ch}.filter(ref,data_bckup_raw,0.05^nan_loop_ctr*StepSize_raw,a_logical,reset_flag); % start over from backup copy
            reset_flag = false;            
            [~, ~, ~] = adaptive_filt_cell_raw{i_ch,j_ch}.filter(ref,data_bckup_raw,0.05^nan_loop_ctr,a_logical,reset_flag); % start over from backup copy
            [~, data_raw, ~] = adaptive_filt_cell_raw{i_ch,j_ch}.filter(ref,data_bckup_raw,0.05^nan_loop_ctr,a_logical,reset_flag); % start over from backup copy
            
            if nan_loop_ctr > 10
                error('cant fix nans')
            end
            
        end 
        
        data_bckup_raw = data_raw; % it made it through, so update the backup
        
    end
    
    if n_read == 2 % bp filter raw data after the artifact removal
        [data_raw,zf0_raw] = filter(b_H_raw,a_H_raw,data_raw);
        [data_raw,zf01_raw] = filter(b_S,a_S,data_raw);
        [data_raw,zf02_raw] = filter(b_L_raw,a_L_raw,data_raw);
    else
        [data_raw,zf0_raw] = filter(b_H_raw,a_H_raw,data_raw,zf0_raw);
        [data_raw,zf01_raw] = filter(b_S,a_S,data_raw,zf01_raw);
        [data_raw,zf02_raw] = filter(b_L_raw,a_L_raw,data_raw,zf02_raw);
    end
   
%     unpadded_data = data(prev_data_chunk_len+1:end,:);
%     art_rem_data(:,i_ch) = unpadded_data;  % should add a lpf before decimation
%     unpadded_data = data(prev_data_chunk_len+1:end,:);
    art_rem_data(:,i_ch) = data(end-size(art_rem_data,1)+1:end,:);  % should add a lpf before decimation
    art_rem_data_raw(:,i_ch) = data_raw(end-size(art_rem_data_raw,1)+1:end,:);  % should add a lpf before decimation
end

previous_data_chunk = event.Data(end-prev_data_chunk_len:end,:); % save small chunk of the end of the data for next time to help initilizae adaptive filters
if any(isnan(previous_data_chunk))
    warning('event Data contained NaNs!')
    previous_data_chunk = zeros(size(previous_data_chunk));
end

%%
if any(isnan(art_rem_data))
    error('lms returning nans')
end
if any(isnan(art_rem_data_raw))
    nan_chs = num2str(find(sum(isnan(at_rem_data_raw),1)));%determine which channel is producing NaNs
    error(['lms raw returning nan, ch = ' nan_chs])
end

%% save art rem data
%EFF21 - no longer check size from matfile but keep track in few vars
[r_ard, c_ard] = size(art_rem_data);
[r_ard_raw, c_ard_raw] = size(art_rem_data_raw);
if DoSaveData    
    if co_ard~=c_ard
        error('data dims not equal in columns')
    end

    if co_ard~=c_ard_raw
        error('data dims not equal in columns')
    end

    mf2.art_rem_data(ro_ard+(1:r_ard), :) =single(art_rem_data); %EFF21, remove time and turn to single
    mf2.art_rem_data_raw(ro_ard+(1:r_ard_raw), :) =single(art_rem_data_raw); %EFF21, remove time and turn to single
    ro_ard=ro_ard+r_ard; %EFF21? Assumed that art rem data and art rem data raw always have same number of rows
end

deci_freq = 5000; %%%CKM??? hardcoded deci_freq=5000, why 5000

% [b_deci a_deci] = butter(3,0.425*deci_freq/(fs_rec/2),'low');
% global zfn1
% if n_read == 2
%     [art_rem_data, zfn1] = filter(b_deci,a_deci,art_rem_data);
% else
%     [art_rem_data, zfn1] = filter(b_deci,a_deci,art_rem_data,zfn1);
% end


%% plot decimated data
deci = fix(fs_rec/deci_freq); % downsample plots to about 5 khz
fs_rec_deci = fs_rec/deci;

% MainPlot_handle = findobj('Tag', 'A_MainPlot');

data_deci = event.Data(1:deci:end,:);

data_deci(:,detect_ch_index) = art_rem_data(1:deci:end,:);
art_rem_data_raw = art_rem_data_raw(1:deci:end,:);
time_deci = event.TimeStamps(1:deci:end);
n_t_deci = length(time_deci);

%load in settings to position the data vertically within the plot
channel_scaling = str2num(get(handles.ET_channel_spacing,'String'));
channel_spacing = [1:n_ch-n_ch_out 2*(1:n_ch_out)];

%refresh the plot to avoid loading in too much data slowing things down for
%real time detection
t1 = round(time_deci(1));
time_x = round(str2num(get(handles.ET_time_plot,'String')));    %#s to display on x axis of plot
if ~DoExisting
    if mod(t1,time_x) == 0                      %until window is filled with data, just add graph to existing window
        hold(plot_handle, 'off')
    end
else
    if get(handles.CheckFluid,'value')
      if mod(t1,time_x) == 0                      %until window is filled with data, just add graph to existing window
        hold(plot_handle, 'off')
      end 
    end
end

ax = findall(plot_handle,'type', 'axes');   %get the axes object 
ax.ColorOrderIndex = 1;                     %Make sure colors of graphs are consistent by resetting order

%% plot raw data
raw_deci = 1;
raw_d = event.Data(1:raw_deci:end,:);   %take raw data (multiple columns)
raw_d = raw_d - mean(raw_d,1);          %remove mean for each column
plot(plot_handle, event.TimeStamps(1:raw_deci:end), raw_d*diag(1./channel_scaling)+repmat(channel_spacing, length(event.TimeStamps(1:raw_deci:end)),1),'LineWidth',1)
hold(plot_handle, 'on')
%%

% commented out for video
% plot(plot_handle, time_deci, data_deci*diag(1./channel_scaling)+repmat(channel_spacing, n_t_deci,1),'LineWidth',1)
% hold(plot_handle, 'on')
% plot(plot_handle, event.TimeStamps, art_rem_data*diag(1./channel_scaling((1:n_ch_out)*2-1))+repmat(channel_spacing((1:n_ch_out)*2-1), length(event.TimeStamps),1),'k')

%Set the limits of the plot
xlim(plot_handle, [t1-rem(t1,time_x) t1-rem(t1,time_x)+time_x])
ylim(plot_handle,[-1 max(channel_spacing(:))+2])


%% seizure detection code
% decimate and remove mean

global prev_digi_data prev_digi_time
digi_data = data_deci(:,digital_ch_index);
digi_time = time_deci(:);

digi_data_concat = [prev_digi_data; digi_data]; %prev_digi_data is initialized at 100 zeros
digi_time_concat = [prev_digi_time; digi_time];

if length(digital_ch_index) ~= n_ch_out
    error('n_ch_out and number of digital channels must be equal')
end

%% ================================
hpf = 500;
lpf = 1500;

[b_hpf a_hpf] = butter(3,hpf/(fs_rec_deci/2),'high');
[b_lpf a_lpf] = fir1(3,lpf/(fs_rec_deci/2),'low');

[b_hpf_raw a_hpf_raw] = butter(3,1/(fs_rec_deci/2),'high');
[b_lpf_raw a_lpf_raw] = butter(3,150/(fs_rec_deci/2),'low');

% % [b_hpf a_hpf] = fir1(20,hpf/(fs_rec_deci/2),'high');
% % [b_lpf a_lpf] = fir1(7,lpf/(fs_rec_deci/2),'low');
% % 
% % [b_hpf_raw a_hpf_raw] = butter(3,1/(fs_rec_deci/2),'high');
% % [b_lpf_raw a_lpf_raw] = butter(3,150/(fs_rec_deci/2),'low');

% data_deci = data_deci - repmat(mean(data_deci,1),size(data_deci,1),1); % remove mean
% can't remove mean of digital channels, hpf should handle this

global zf0p5 zf1 zf2 zf3 zf1_raw zf2_raw

if n_read == 2 %n_read=2 on the first iteration of emptying the buffer
%     [filt_data zf1] = filter(b_lpf,1,filter(b_hpf,1,data_deci(:,detect_ch_logical)));
    [filt_data zf0p5] = filter(b_hpf,a_hpf,data_deci(:,detect_ch_logical));
    [filt_data zf1] = filter(b_lpf,a_lpf,filt_data);
    filt_data= abs(filt_data);
    envf = 25;
    [b_env, a_env] = butter(3,envf/(fs_rec_deci/2),'low');
    [filt_data zf2] = filter(b_env,a_env,filt_data);
    envf_hp = 2;
    [b_env_hp a_env_hp ]= butter(3,envf_hp/(fs_rec_deci/2),'high');
    [filt_data zf3] = filter(b_env_hp,a_env_hp,filt_data);
    
    [art_rem_data_raw, zf1_raw] = filter(b_hpf_raw, a_hpf_raw, art_rem_data_raw);
    [art_rem_data_raw, zf2_raw] = filter(b_lpf_raw, a_lpf_raw, art_rem_data_raw);
    
else % if n_read>2, then can use filter delays coefs (starting at second iteration
    [filt_data zf0p5] = filter(b_hpf,a_hpf,data_deci(:,detect_ch_logical),zf0p5);
    [filt_data zf1] = filter(b_lpf,a_lpf,filt_data,zf1);
    filt_data= abs(filt_data);
    envf = 25;
    [b_env, a_env] = butter(3,envf/(fs_rec_deci/2),'low');
    [filt_data zf2] = filter(b_env,a_env,filt_data,zf2);
    envf_hp = 2;
    [b_env_hp a_env_hp ]= butter(3,envf_hp/(fs_rec_deci/2),'high');
    [filt_data zf3] = filter(b_env_hp,a_env_hp,filt_data,zf3);
    
    [art_rem_data_raw, zf1_raw] = filter(b_hpf_raw, a_hpf_raw, art_rem_data_raw, zf1_raw);
    [art_rem_data_raw, zf2_raw] = filter(b_lpf_raw, a_lpf_raw, art_rem_data_raw, zf2_raw);
end
%%%%%%%%%%%%%

%%

detect_deci = fix(fs_rec_deci/5000); %hardcoded to be the same as deci

%EFF21 TODO - can move this out of loop if it adds up
if n_read == 2 && DoSaveToFile %first iteration of buffer data (seconds)
   mf.detect_deci = detect_deci; 
end

fs_deci = fs_rec_deci/detect_deci;
dt_deci = 1/fs_deci;
%note that both art_rem_data and art_rem_data_raw are DIFFERENT then when
%they were saved to file. So detect_data is not simply reassignment that
%could be skipped.
detect_data = filt_data(1:detect_deci:end,:);    % decimate
detect_data = detect_data - mean(detect_data,1); % remove mean

detect_data_raw = art_rem_data_raw(1:detect_deci:end,:);     % decimate
detect_data_raw = detect_data_raw - mean(detect_data_raw,1); % remove mean

n_t_detect_deci = size(detect_data,1);
detect_time_deci = time_deci(1:detect_deci:end);

plot(plot_handle, detect_time_deci, 50*detect_data*diag(1./channel_scaling(detect_ch_logical))+repmat(channel_spacing(detect_ch_logical), n_t_detect_deci,1))
plot(plot_handle, detect_time_deci, 1*detect_data_raw*diag(1./channel_scaling(detect_ch_logical))+repmat(channel_spacing(detect_ch_logical), n_t_detect_deci,1),'Color',[0.5 0.5 0.5])



%% save decimated and filtered data
[r, ~] = size(detect_data);
%[r_raw, c_raw] = size(detect_data_raw);
if DoSaveData
    %EFF21 - no longer check size from file, but keep track with ro_dd and co_dd
    mf2.detect_data(ro_dd+(1:r), :) = single(detect_data);         %EFF21 (remove time)(converto to single)
    mf2.detect_data_raw(ro_dd+(1:r), :) = single(detect_data_raw); %EFF21 (remove time)(converto to single)
    ro_dd=ro_dd+r;
end

%% spike analysis
dist_min_vec = str2num(get(handles.ET_spike_dist_min,'String')); % ms?
dist_max_vec = str2num(get(handles.ET_spike_dist_max,'String')); % ms?
SpikesUpper_vec = str2num(get(handles.ET_threshold_pos,'String')); % SD?
SpikesLower_vec = str2num(get(handles.ET_threshold_neg,'String')); % SD?
SpikesPosMax_vec = str2num(get(handles.ET_maxAmp_pos,'String')); % SD?
SpikesNegMax_vec = str2num(get(handles.ET_maxAmp_neg,'String')); % SD?
ht_perc_vec = str2num(get(handles.ET_width_at_percent_height,'String')); %Width of a spike is determined at specified height of spikes (fraction)
req_width_vec = str2num(get(handles.ET_min_width,'String')); %Spike is considered wide if width exceeds this threshold (in s)
req_width_vec_max = str2num(get(handles.ET_max_width,'String')); %Spike is considered too wide if width exceeds this threshold (in s)

dist_min_vec_raw = str2num(get(handles.ET_spike_dist_min_raw,'String')); % ms?
dist_max_vec_raw = str2num(get(handles.ET_spike_dist_max_raw,'String')); % ms?
SpikesUpper_vec_raw = str2num(get(handles.ET_threshold_pos_raw,'String')); % SD?
SpikesLower_vec_raw = str2num(get(handles.ET_threshold_neg_raw,'String'));% SD?
SpikesPosMax_vec_raw = str2num(get(handles.ET_maxAmp_raw_pos,'String')); % SD?
SpikesNegMax_vec_raw = str2num(get(handles.ET_maxAmp_raw_neg,'String')); % SD?
ht_perc_vec_raw = str2num(get(handles.ET_width_at_percent_height_raw,'String')); %Width of a spike is determined at specified height of spikes (fraction)
req_width_vec_raw = str2num(get(handles.ET_min_width_raw,'String')); %Spike is considered wide if width exceeds this threshold (in s)
req_width_vec_raw_max = str2num(get(handles.ET_max_width_raw,'String')); %Spike is considered too wide if width exceeds this threshold (in s)

%EFF21 changes to spike count history - removed 3d structure
spike_count_history = circshift(spike_count_history,-1); % shift values up
spike_count_history_raw = circshift(spike_count_history_raw,-1); % shift values up
spike_count_history_combined = circshift(spike_count_history_combined,-1); % shift values up
spike_times_buffer_cell = circshift(spike_times_buffer_cell,-1);
spikes_trigger=circshift(spikes_trigger,-1);
spikes_trigger_stop=circshift(spikes_trigger_stop,-1);

for i_ch = 1:n_ch_out
    Settings.dist_min = dist_min_vec(i_ch);
    Settings.dist_max = dist_max_vec(i_ch);
    Settings.SpikesUpper = SpikesUpper_vec(i_ch);
    Settings.SpikesLower = SpikesLower_vec(i_ch);
    Settings.SpikesPosMax = SpikesPosMax_vec(i_ch);
    Settings.SpikesNegMax = SpikesNegMax_vec(i_ch);
    Settings.ht_perc = ht_perc_vec(i_ch);
    Settings.req_width = req_width_vec(i_ch);
    Settings.req_width_max = req_width_vec_max(i_ch);
    
    Settings_raw.dist_min = dist_min_vec_raw(i_ch);
    Settings_raw.dist_max = dist_max_vec_raw(i_ch);
    Settings_raw.SpikesUpper = SpikesUpper_vec_raw(i_ch);
    Settings_raw.SpikesLower = SpikesLower_vec_raw(i_ch);
    Settings_raw.SpikesPosMax = SpikesPosMax_vec_raw(i_ch);
    Settings_raw.SpikesNegMax = SpikesNegMax_vec_raw(i_ch);
    Settings_raw.ht_perc = ht_perc_vec_raw(i_ch);
    Settings_raw.req_width = req_width_vec_raw(i_ch);
    Settings_raw.req_width_max = req_width_vec_raw_max(i_ch);
    
    %look at the called function to see which other vars were calculated.
    %Adjusted functions to only calculate and return vars we wanted. For
    %future changes find original version of the functions
    [posspikes_wide,negspikes_wide]=EffSpikeFinder(detect_data(:,i_ch),fs_deci,Settings);
    [posspikes_wide_raw,negspikes_wide_raw]=EffSpikeFinder(detect_data_raw(:,i_ch),fs_deci,Settings_raw);

    n_spikes_this_chunk = size(negspikes_wide,1)+size(posspikes_wide,1);
    n_spikes_this_chunk_raw = size(negspikes_wide_raw,1)+size(posspikes_wide_raw,1);
    
    spike_count_history(10,i_ch) = n_spikes_this_chunk; %HARDCODED
    spike_count_history_raw(10,i_ch) = n_spikes_this_chunk_raw; %HARDCODED
    
    stim_start_concat_ind = find(diff(digi_data_concat(:,i_ch))>0.5);
    stim_stop_concat_ind = find(diff(digi_data_concat(:,i_ch))<-0.5);
    stim_start_time = digi_time_concat(stim_start_concat_ind);
    stim_stop_time = digi_time_concat(stim_stop_concat_ind);
    
    blank_time = 0.2; % hard coded one fifth of a second
    
    %Process wide positive spikes
    if not(isempty(posspikes_wide))
        pos_spike_times = detect_time_deci(posspikes_wide(:,1));
        pos_spike_val = posspikes_wide(:,2);
        
        % check for and remove spikes immediately after the start and stop of the stimulus
        start_spikes_to_blank_logical = false(size(pos_spike_times));
        stop_spikes_to_blank_logical = false(size(pos_spike_times));
        combined_spikes_to_blank_logical = false(size(pos_spike_times));
        
        if not(isempty(stim_start_time))
            start_spikes_to_blank_logical = and( stim_start_time<pos_spike_times, pos_spike_times < stim_start_time+blank_time );
        end
      
        combined_spikes_to_blank_logical = or(start_spikes_to_blank_logical, stop_spikes_to_blank_logical); % this was the bug, should be OR, not AND!!!!
        
        pos_spike_times(combined_spikes_to_blank_logical) = [];
        pos_spike_val(combined_spikes_to_blank_logical) = [];
 
        n_pos_spike(1,i_ch) = length(pos_spike_times);
        
        if not(isempty(pos_spike_times)) && DoSaveToFile
            %EFF21 - store times and val in same array to reduce #writes to file
            mf.pos_spikes(pos_spike_count(1,i_ch)+uint32((1:n_pos_spike(1,i_ch))),2*i_ch-1:2*i_ch) = [pos_spike_times pos_spike_val];
        end
        pos_spike_count(1,i_ch) = pos_spike_count(1,i_ch) + uint32(n_pos_spike(1,i_ch));
    else
        pos_spike_times = [];
    end

    %Process wide negative spikes
    if not(isempty(negspikes_wide))
        neg_spike_times = detect_time_deci(negspikes_wide(:,1));
        neg_spike_val = negspikes_wide(:,2);
         
        % check for and remove spikes immediately after the start and stop of the stimulus
        start_spikes_to_blank_logical = false(size(neg_spike_times));
        stop_spikes_to_blank_logical = false(size(neg_spike_times));
        combined_spikes_to_blank_logical = false(size(neg_spike_times));
        
        if not(isempty(stim_start_time))
            start_spikes_to_blank_logical = and( stim_start_time<neg_spike_times, neg_spike_times < stim_start_time+blank_time );
        end
        
        if not(isempty(stim_stop_time))
            stop_spikes_to_blank_logical = and( stim_stop_time<neg_spike_times, neg_spike_times < stim_stop_time+blank_time );
        end
        
        combined_spikes_to_blank_logical = or(start_spikes_to_blank_logical, stop_spikes_to_blank_logical);

        neg_spike_times(combined_spikes_to_blank_logical) = [];
        neg_spike_val(combined_spikes_to_blank_logical) = [];
        
        n_neg_spike(1,i_ch) = length(neg_spike_times); 
        if DoSaveToFile
            %EFF21 - store times and val in same array to reduce #writes to file
            mf.neg_spikes(neg_spike_count(1,i_ch)+uint32((1:n_neg_spike(1,i_ch))),2*i_ch-1:2*i_ch) = [neg_spike_times neg_spike_val];
        end
        neg_spike_count(1,i_ch) = neg_spike_count(1,i_ch) + uint32(n_neg_spike(1,i_ch));
    else
        neg_spike_times = [];
    end
    
    %% raw spikes
    
    %deal with positive wide spikes in raw detection
    if not(isempty(posspikes_wide_raw))
        pos_spike_times_raw = detect_time_deci(posspikes_wide_raw(:,1));
        pos_spike_val_raw = posspikes_wide_raw(:,2);
        
        % check for and remove spikes immediately after the start and stop of the stimulus
        start_spikes_to_blank_logical = false(size(pos_spike_times_raw));
        stop_spikes_to_blank_logical = false(size(pos_spike_times_raw));
        combined_spikes_to_blank_logical = false(size(pos_spike_times_raw));
        
        if not(isempty(stim_start_time))
            start_spikes_to_blank_logical = and( stim_start_time<pos_spike_times_raw, pos_spike_times_raw < stim_start_time+blank_time );
        end
        
        
        combined_spikes_to_blank_logical = or(start_spikes_to_blank_logical, stop_spikes_to_blank_logical); % this was the bug, should be OR, not AND!!!!
        
        pos_spike_times_raw(combined_spikes_to_blank_logical) = [];
        pos_spike_val_raw(combined_spikes_to_blank_logical) = [];
 
        n_pos_spike_raw(1,i_ch) = length(pos_spike_times_raw);
        
        if not(isempty(pos_spike_times_raw)) && DoSaveToFile
            %EFF21 - store times and val in same array to reduce #writes to file
            mf.pos_spikes_raw(pos_spike_count_raw(1,i_ch)+uint32((1:n_pos_spike_raw(1,i_ch))),2*i_ch-1:2*i_ch) = [pos_spike_times_raw pos_spike_val_raw];
        end
        pos_spike_count_raw(1,i_ch) = pos_spike_count_raw(1,i_ch) + uint32(n_pos_spike_raw(1,i_ch));
    else
        pos_spike_times_raw = [];
    end

    %deal with negative wide spikes in raw detection
    if not(isempty(negspikes_wide_raw))
        neg_spike_times_raw = detect_time_deci(negspikes_wide_raw(:,1));
        neg_spike_val_raw = negspikes_wide_raw(:,2);
         
        % check for and remove spikes immediately after the start and stop of the stimulus
        start_spikes_to_blank_logical = false(size(neg_spike_times_raw));
        stop_spikes_to_blank_logical = false(size(neg_spike_times_raw));
        combined_spikes_to_blank_logical = false(size(neg_spike_times_raw));
        
        if not(isempty(stim_start_time))
            start_spikes_to_blank_logical = and( stim_start_time<neg_spike_times_raw, neg_spike_times_raw < stim_start_time+blank_time );
        end
        
        if not(isempty(stim_stop_time))
            stop_spikes_to_blank_logical = and( stim_stop_time<neg_spike_times_raw, neg_spike_times_raw < stim_stop_time+blank_time );
        end
        
        combined_spikes_to_blank_logical = or(start_spikes_to_blank_logical, stop_spikes_to_blank_logical);

%         if any(combined_spikes_to_blank_logical)
%             display('stim spike blanked!')
%         end
        
        neg_spike_times_raw(combined_spikes_to_blank_logical) = [];
        neg_spike_val_raw(combined_spikes_to_blank_logical) = [];
        
        n_neg_spike_raw(1,i_ch) = length(neg_spike_times_raw); 
        if not(isempty(neg_spike_times_raw))&& DoSaveToFile
            %EFF21 - store times and val in same array to reduce #writes to file
            mf.neg_spikes_raw(neg_spike_count_raw(1,i_ch)+uint32((1:n_neg_spike_raw(1,i_ch))),2*i_ch-1:2*i_ch) = [neg_spike_times_raw neg_spike_val_raw];
        end
        neg_spike_count_raw(1,i_ch) = neg_spike_count_raw(1,i_ch) + uint32(n_neg_spike_raw(1,i_ch));
    else
        neg_spike_times_raw = [];
    end
    
    
    %% spike logic
    if i_ch == 1
        logic_vec = str2num(get(handles.ET_spike_logic_1,'String'));
    elseif i_ch == 2
        logic_vec = str2num(get(handles.ET_spike_logic_2,'String'));
    elseif i_ch == 3
        logic_vec = str2num(get(handles.ET_spike_logic_3,'String'));
    elseif i_ch == 4
        logic_vec = str2num(get(handles.ET_spike_logic_4,'String'));
    end
    
    [~, cross, logic_val] = find(logic_vec);
    if length(cross)~= 1
        error('only one logic statement per mouse, others must be zero')
    end
    
    if cross == 1 % p x n
        times_A = pos_spike_times(:);
        times_B = neg_spike_times(:);
    elseif cross == 2 % p x p_raw
        times_A = pos_spike_times(:);
        times_B = pos_spike_times_raw(:);
    elseif cross == 3 % p x n_raw
        times_A = pos_spike_times(:);
        times_B = neg_spike_times_raw(:);
    elseif cross == 4 % n x p_raw
        times_A = neg_spike_times(:);
        times_B = pos_spike_times_raw(:);
    elseif cross == 5 % n x n_raw
        times_A = neg_spike_times(:);
        times_B = neg_spike_times_raw(:);
    elseif cross == 6 % p_raw x n_raw
        times_A = pos_spike_times_raw(:);
        times_B = neg_spike_times_raw(:);
    elseif cross == 7 % p x p_raw or n_raw
        times_A = pos_spike_times(:);
        times_B = pos_spike_times_raw(:);
        times_C = neg_spike_times_raw(:);
    else
        error('value outside 1:7, one logic value of 1 to 7 must be set per mouse')
    end
    
    coincident_spike_times = [];
    
    if and(cross == 7, logic_val ~= 5)
        error('if using position 7, it must be logic value 5 for and(p_hf or(p_raw, n_raw))')
    end
    
    if logic_val == 1 % or
        
        if isempty(times_A)
            times_A = NaN;
        end
             
        if isempty(times_B)
            times_B = NaN;
        end
        
        and_window = str2num(get(handles.ET_and_window,'String'));
                
        coincident_spike_logical = min(pdist2(times_B,times_A,'cityblock'),[],2)<and_window(i_ch);
        times_B(coincident_spike_logical) = []; % remove any Bs that are coincident with A to avoid double counts
        
        coincident_spike_times = [times_A; times_B];
        coincident_spike_times(isnan(coincident_spike_times)) = []; % remove NaNs
        
        valid_spike_times = coincident_spike_times;
               
    elseif logic_val == 2 % and
        
        and_window = str2num(get(handles.ET_and_window,'String'));
        
        if isempty(times_A)
            times_A = NaN;
        end
             
        if isempty(times_B)
            times_B = NaN;
        end
        
        coincident_spike_logical = min(pdist2(times_A,times_B,'cityblock'),[],2)<and_window(i_ch);
        coincident_spike_times = times_A(coincident_spike_logical); 
        coincident_spike_times(isnan(coincident_spike_times)) = [];
        
        valid_spike_times = coincident_spike_times(:);
        
    elseif logic_val == 3 % and not
        
        and_window = str2num(get(handles.ET_and_window,'String'));
        
        times_A_cpy = times_A;
        
        if isempty(times_A)
            times_A = NaN;
        end
             
        if isempty(times_B)
            times_B = NaN;
        end
        
%         if not(isempty(times_A_cpy))
%             error('stop!')
%         end
        
        coincident_spike_logical = min(pdist2(times_A,times_B,'cityblock'),[],2)<and_window(i_ch);
        times_A(coincident_spike_logical) = []; 
        times_A(isnan(times_A)) = [];
        
        valid_spike_times =  times_A(:);
        
    elseif logic_val == 4 % not and
        
        and_window = str2num(get(handles.ET_and_window,'String'));
        
        if isempty(times_A)
            times_A = NaN;
        end
             
        if isempty(times_B)
            times_B = NaN;
        end
        
        coincident_spike_logical = min(pdist2(times_B,times_A,'cityblock'),[],2)<and_window(i_ch);
        times_B(coincident_spike_logical) = [];  
        times_B(isnan(times_B)) = [];
        
        valid_spike_times = times_B(:);
        
    elseif logic_val == 5 % and(p_hf, or(p_raw,n_raw))
        
        and_window = str2num(get(handles.ET_and_window,'String'));
        
        if isempty(times_A)
            times_A = NaN;
        end    
        
        times_B = [times_B; times_C]; % this makes the or() work
        
        if isempty(times_B)
            times_B = NaN;
        end

        coincident_spike_logical = min(pdist2(times_A,times_B,'cityblock'),[],2)<and_window(i_ch);
        coincident_spike_times = times_A(coincident_spike_logical); 
        coincident_spike_times(isnan(coincident_spike_times)) = [];
        
        valid_spike_times = coincident_spike_times(:);
        
    else
        error('logic value must be in range 1 to 4')
    end
    
    valid_spike_times = sort(valid_spike_times);
    if not(isempty(valid_spike_times))
        n_valid_spike(1,i_ch) = length(valid_spike_times);
        valid_spike_count(1,i_ch) = valid_spike_count(1,i_ch) + n_valid_spike(1,i_ch);
        if DoSaveToFile
            mf.valid_spike_times(valid_spike_count(1,i_ch)+uint32((1:n_valid_spike(1,i_ch))),i_ch) = valid_spike_times;
        end
    end
    
    spike_times_buffer_cell{end,i_ch} = valid_spike_times(:);
    spike_times_buffer_timeCombined{1,i_ch} = sort( cell2mat( spike_times_buffer_cell(:,i_ch) ) );
    spike_count_history_combined(10,i_ch) = length(valid_spike_times); %HARDCODED
    
    %%
    ch_ind = find(seizure_detection_ch_vec(i_ch) == ch_in_vec);
    if ~isempty(pos_spike_times)    
        plot(plot_handle, pos_spike_times,(1./channel_scaling(ch_ind))*50*pos_spike_val+ch_ind,'*r')
    end
    
    if ~isempty(neg_spike_times)
        plot(plot_handle, neg_spike_times,(1./channel_scaling(ch_ind))*50*neg_spike_val+ch_ind,'*r')
    end
    
    if ~isempty(pos_spike_times_raw)
        plot(plot_handle, pos_spike_times_raw,(1./channel_scaling(ch_ind))*pos_spike_val_raw+ch_ind,'*g')
    end
    
    if ~isempty(neg_spike_times_raw)
        plot(plot_handle, neg_spike_times_raw,(1./channel_scaling(ch_ind))*neg_spike_val_raw+ch_ind,'*b')
    end
    
    if ~isempty(valid_spike_times)
        plot(plot_handle, valid_spike_times,ch_ind,'*k')
    end
end

min_spikes_per_two_s_vec = str2num(get(handles.ET_min_spikes_per_two_s,'String'));
min_spikes_per_two_s_vec_stop = str2num(get(handles.ET_min_spikes_per_two_s_stop,'String'));

start_seizure_n_seconds = 2; % seconds, hard coded

for i_ch = 1:n_ch_out
    s_spikes = sum(spike_count_history_combined(10-start_seizure_n_seconds+1:10,i_ch)); % number of spikes in the last 2 seconds %%%???CKM suppressed output %HARDCODED10
    
    if s_spikes>=min_spikes_per_two_s_vec(i_ch)
        spikes_trigger(10,i_ch) = true;
        Tr_On_Index=[Tr_On_Index; n_read i_ch]; %EFF21 - prepping for sparse indexing
    else           
        spikes_trigger(10,i_ch) = false;
    end
    
    if s_spikes>=min_spikes_per_two_s_vec_stop(i_ch)
        spikes_trigger_stop(10,i_ch) = true;
        Tr_Off_Index=[Tr_Off_Index; n_read i_ch]; %EFF21 - prepping for sparse indexing
    else           
        spikes_trigger_stop(10,i_ch) = false;
    end
end
%EFF21 - can get rid of else part if initialized properly.
if DoSaveToFile
    %EFF21: CHANGED STRUCTURE to 2d array (used to be that spike_count_history(:,:,2) is
    %spike_count_history(:,:,2) with everything shifted down 1).
    %MAYBE save every 10 to reduce save time, but then must add in final
    %save as well.
    if mod(n_read,10)==0  %10 IS HARDCODED: CHANGE IF BLOCK IS INITIALIZED DIFFERENTLY
        mf.spike_count(n_read-9:n_read,:) = [spike_count_history ...        
                                                     spike_count_history_raw ...
                                                     spike_count_history_combined];
    end
    
    %EFF21 TODO: don't save each time, but rather save at certain
    %intervals, using sparse structure HARDCODED
    %mf.spikes_trigger(n_read,:) = spikes_trigger(10,:);
    %mf.spikes_trigger_stop(n_read,:) = spikes_trigger_stop(10,:);
end

% for i_ch = 1:n_ch_out
%     if spikes_trigger(n_read,i_ch) == 1
%         spike_times_buffer_cell{end,i_ch};
%         if not(isempty(spike_times_buffer_cell{end,i_ch}))
%             last_spike_time(1,i_ch) = max(spike_times_buffer_cell{end,i_ch});
%         end
%         if size(last_spike_time,1)>1
%             error('spikes_triggered not 1x2')
%         end
%     end
% end

for i_ch = 1:n_ch_out
    if spikes_trigger_stop(10,i_ch) == 1 %EFF21 TODO - change spikes_trigger_stop and spikes_trigger to small arrays that get looped
        %spike_times_buffer_cell{end,i_ch}; %%EFF21 NOT DOING ANYTHING
        if not(isempty(spike_times_buffer_cell{end,i_ch}))
            last_spike_time(1,i_ch) = max(spike_times_buffer_cell{end,i_ch});
        end
        if size(last_spike_time,1)>1
            error('spikes_triggered not 1xn')
        end
    end
end


%% seizure starts
%ready_for_next_sz = ready_for_next_sz %EFF21 CKM??? What is the point of this line, besides communication
disp(['ready_for_next_sz = ' num2str(ready_for_next_sz)]);
%EFF21, reduced role of seizure_On, only storing LAST ENTRY, NEVER saving
Seizure_On(2,:)=Seizure_On(1,:);   %row2: last index     %EFF21: getting ready for next iteration
                                   %row1: current index, initialize to previous index, then will see what else changes in next steps
if n_read > 5
    for i_ch = 1:n_ch_out
        if ((time_out(1,i_ch) <=0) && (ready_for_next_sz(i_ch) == true) && (spikes_trigger(10,i_ch) == 1) && (spikes_trigger(9,i_ch) == 1) && (spikes_trigger(8,i_ch) == 1) && (spikes_trigger(7,i_ch) == 0) && (spikes_trigger(6,i_ch) == 0) && (spikes_trigger(5,i_ch) == 0) && (spikes_trigger(4,i_ch) == 0) && (spikes_trigger(3,i_ch) == 0) && (spikes_trigger(2,i_ch) == 0)) % if seizure conditions satisfied, and time_out has counted down  %HARDCODED EFF21
            Seizure_On(1,i_ch) = 1; %seizure is happening currently
        %else % no change
            %Seizure_On(1,i_ch) = Seizure_On(1,i_ch);
        end
    %else
    %    Seizure_On(1,i_ch) = 0;
    end
end


%% seizure ends
if n_read>5
    for i_ch = 1:n_ch_out
%         if spikes_trigger(n_read,i_ch) == 0 && spikes_trigger(n_read-1,i_ch) == 0 && spikes_trigger(n_read-2,i_ch) == 0 && spikes_trigger(n_read-3,i_ch) == 0 && spikes_trigger(n_read-4,i_ch) == 0
        if spikes_trigger_stop(10,i_ch) == 0 && spikes_trigger_stop(9,i_ch) == 0 && spikes_trigger_stop(8,i_ch) == 0 % && spikes_trigger(n_read-3,i_ch) == 0 && spikes_trigger(n_read-4,i_ch) == 0
            Seizure_On(1,i_ch) = 0; %no seizure currently
        end
    end
%else %EFF21: not necessary since it is initialized to false
%    Seizure_On(1,:) = 0; %no seizure currently
end

%if n_read > 1 %EFF21 - not necessary since we don't index to negs anymore
Seizure_Off = and(Seizure_On(1,:)==0, Seizure_On(2,:)==1); %turns off if there was seizure, but not anymore
%else
%    Seizure_Off(n_read,:) = 0;
%end

%% stim when seizure detected
time_out = time_out - in_chunk;                             %time_out counts down from ET_excusion_time.  If it is at or below zero, then you can stim.
stim_flag = and(Seizure_On(1,:)==1, Seizure_On(2,:)==0);    %turns on when seizure starts
Seizure_Start_Ind(stim_flag) = n_read;
Seizure_Count = Seizure_Count + stim_flag;                  %keep track of seizure numbers across all files
Seizure_Count_Current = Seizure_Count_Current + stim_flag;  %keep track of seizure numbers in current file


%% log Seizure_On and Seizure_Off
if DoSaveToFile
    %mf.Seizure_On(n_read,:) = Seizure_On(n_read,:); %EFF21: REMOVE VAR,UNNECESSARY
    %mf.Seizure_Off(n_read,:) = Seizure_Off;    %EFF21 - better to storein mem and save at end
    %mf.stim_flag(n_read,:) = stim_flag;
    %mf.time_out(n_read,:) = time_out;          %EFF21 - no need to save - can be obtained from Seizure_On, Seizure_Off and settings on how to initialize
    
    %Save index of the stimulation starts (will be saved to file later)
    szlocs=find(stim_flag(1,:));
    for i_l=1:length(szlocs)
        Sz_On_Index=[Sz_On_Index; n_read szlocs(i_l)];
    end
    
    %Save index of the stimulation ends (will be saved to file later)
    szlocs=find(Seizure_Off(1,:));
    for i_l=1:length(szlocs)
        Sz_Off_Index=[Sz_Off_Index; n_read szlocs(i_l)];
    end
    
end

for i_ch = 1:n_ch_out
    if stim_flag(i_ch) == 1 %stim flag has (1s) chunks that are flagged for stim
        Seizure_Start_First_Spike_Time(1,i_ch) = min(spike_times_buffer_timeCombined{1,i_ch});
        Seizure_Start_First_Spike_Time_Current(1,i_ch) = Seizure_Start_First_Spike_Time(1,i_ch);
        if DoExisting
            plot(plot_handle, Seizure_Start_First_Spike_Time(1, i_ch) ,(channel_spacing(2*i_ch-1)),'>y');
        end
        if DoSaveToFile
            mf.Seizure_Start_First_Spike_Time(Seizure_Count(i_ch), i_ch) = Seizure_Start_First_Spike_Time(1, i_ch);
            if DoSaveData
                mf2.Seizure_Start_First_Spike_Time_Current(Seizure_Count_Current(i_ch), i_ch) = Seizure_Start_First_Spike_Time_Current(1, i_ch);
            end
        end
    end
end

%% count seizures
if DoExisting %Alternative code to use when doing reanalysis of existing data, rather than online analysis
    for i_ch_out=1:n_ch_out
      if Seizure_Off(1,i_ch_out)==1
        time_out(1,i_ch_out) = int32(str2num(get(handles.ET_exclusion_time,'String'))+5*(i_ch_out-1)); % staggers them

        %Determine seizure length
        Seizure_Duration(Seizure_Count(1,i_ch_out),i_ch_out) = last_spike_time(1,i_ch_out) - Seizure_Start_First_Spike_Time(1, i_ch_out);        
        duration_amp_freq(Seizure_Count(1,i_ch_out),(i_ch_out-1)*3+1) = Seizure_Duration(Seizure_Count(1,i_ch_out),i_ch_out);
       
        %If we are in Offline mode (DoExisting =1), then add seizures to
        %graph. plot line Seizure_start_first_spike_time+[0,Seizure_Duration]
        if 1%stim_flag(i_ch_out)
                
            disp(['Ch' num2str(i_ch_out) ': Seizure at: ' num2str(Seizure_Start_First_Spike_Time(1, i_ch_out)) '-' num2str(last_spike_time(1,i_ch_out))]);
            plot(plot_handle, [Seizure_Start_First_Spike_Time(1, i_ch_out) last_spike_time(1,i_ch_out)],(channel_spacing(2*i_ch_out-1))+[0 0],'y');
            plot(plot_handle, last_spike_time(1,i_ch_out) ,(channel_spacing(2*i_ch_out-1)),'<y');
        end    
    
        %Update summary file with info now that seizure has ended
        if DoSaveToFile
            mf.last_spike_time(Seizure_Count(1,i_ch_out),i_ch_out) = last_spike_time(1,i_ch_out);
            mf.Seizure_Duration(Seizure_Count(1,i_ch_out),i_ch_out) = Seizure_Duration(Seizure_Count(1,i_ch_out),i_ch_out);
            mf.duration_amp_freq(Seizure_Count(1,i_ch_out),(i_ch_out-1)*3+1) = Seizure_Duration(Seizure_Count(1,i_ch_out),i_ch_out);
        end
      end
    end
end
%%left off here TODO
if DoSaveData %only do this output part in online mode
for i_ch_out = 1:n_ch_out
    if stim_flag(i_ch_out)
        stim_freq(Seizure_Count(1,i_ch_out),i_ch_out) = next_freq(1,i_ch_out);
        stim_amp(Seizure_Count(1,i_ch_out),i_ch_out) = next_amp(1,i_ch_out);
        stim_width(Seizure_Count(1,i_ch_out),i_ch_out) = next_width(1,i_ch_out);
        stim_ampRatio(Seizure_Count(1,i_ch_out),i_ch_out) = next_ampRatio(1,i_ch_out);
        stim_chgRatio(Seizure_Count(1,i_ch_out),i_ch_out) = next_chgRatio(1,i_ch_out);

        %Update summary file with info now that seizure has started
        %EFF21 TODO: could save this in one save using either structure of
        %many columns IF profile shows that is warranted
        mf.stim_freq(Seizure_Count(1,i_ch_out),i_ch_out) = stim_freq(Seizure_Count(1,i_ch_out),i_ch_out);
        mf.stim_amp(Seizure_Count(1,i_ch_out),i_ch_out) = stim_amp(Seizure_Count(1,i_ch_out),i_ch_out); % mf.stim_amp was not save correctly previously due to typo :( it was ms_ rather than mf., data should be in dur_amp_freq
        mf.stim_width(Seizure_Count(1,i_ch_out),i_ch_out) = stim_width(Seizure_Count(1,i_ch_out),i_ch_out);
        mf.stim_ampRatio(Seizure_Count(1,i_ch_out),i_ch_out) = stim_ampRatio(Seizure_Count(1,i_ch_out),i_ch_out);
        mf.stim_chgRatio(Seizure_Count(1,i_ch_out),i_ch_out) = stim_chgRatio(Seizure_Count(1,i_ch_out),i_ch_out);

        
        duration_amp_freq(Seizure_Count(1,i_ch_out),(i_ch_out-1)*3+2) = stim_amp(Seizure_Count(1,i_ch_out),i_ch_out); % these are redundant and could maybe be removed
        duration_amp_freq(Seizure_Count(1,i_ch_out),(i_ch_out-1)*3+3) = stim_freq(Seizure_Count(1,i_ch_out),i_ch_out); % these are redundant and could maybe be removed
        mf.duration_amp_freq(Seizure_Count(1,i_ch_out),(i_ch_out-1)*3+2) = stim_amp(Seizure_Count(1,i_ch_out),i_ch_out); % these are redundant and could maybe be removed
        mf.duration_amp_freq(Seizure_Count(1,i_ch_out),(i_ch_out-1)*3+3) = stim_freq(Seizure_Count(1,i_ch_out),i_ch_out); % these are redundant and could maybe be removed
        
        %Also save to individual files (CKM) EFF21 - do this in ONE save to save time
        mf2.stim_Current(Seizure_Count_Current(1,i_ch_out),(i_ch_out-1)*5+[1:5]) = [next_freq(1,i_ch_out),...
                                                                                           next_amp(1,i_ch_out),...
                                                                                           next_width(1,i_ch_out),...
                                                                                           next_ampRatio(1,i_ch_out),...
                                                                                           next_chgRatio(1,i_ch_out)];  
    end
end

disp(['Seizure_Count = ' num2str(Seizure_Count)])

freqLo_vec = str2num(get(handles.ET_freqLo,'String'));
freqHi_vec = str2num(get(handles.ET_freqHi,'String'));
ampLo_vec = str2num(get(handles.ET_ampLo,'String'));
ampHi_vec = str2num(get(handles.ET_ampHi,'String'));
widLo_vec = str2num(get(handles.ET_widthLo,'String'));
widHi_vec = str2num(get(handles.ET_widthHi,'String'));
ampRatioLo_vec = str2num(get(handles.ET_ampRatioLo,'String'));
ampRatioHi_vec = str2num(get(handles.ET_ampRatioHi,'String'));
chgRatioLo_vec = str2num(get(handles.ET_chgRatioLo,'String'));
chgRatioHi_vec = str2num(get(handles.ET_chgRatioHi,'String'));

freqBins = str2num(get(handles.ET_freqBins,'String'));
ampBins = str2num(get(handles.ET_ampBins,'String'));
widBins = str2num(get(handles.ET_widthBins,'String'));
ampRatioBins = str2num(get(handles.ET_ampRatioBins,'String'));
chgRatioBins = str2num(get(handles.ET_chgRatioBins,'String'));

optimize_vec = str2num(get(handles.ET_optimize_vec,'String'));

bayes_or_AB_str = get(handles.ET_bayes_or_AB,'String');


%% determine seizure duration
for i_ch_out = 1:n_ch_out
    if Seizure_Off(1,i_ch_out)==1
        time_out(1,i_ch_out) = int32(str2num(get(handles.ET_exclusion_time,'String'))+5*(i_ch_out-1)); % staggers them

        % determine seizure length       
        Seizure_Duration(Seizure_Count(1,i_ch_out),i_ch_out) = last_spike_time(1,i_ch_out) - Seizure_Start_First_Spike_Time(1, i_ch_out);        
        duration_amp_freq(Seizure_Count(1,i_ch_out),(i_ch_out-1)*3+1) = Seizure_Duration(Seizure_Count(1,i_ch_out),i_ch_out);

         mf.last_spike_time(Seizure_Count(1,i_ch_out),i_ch_out) = last_spike_time(1,i_ch_out);
         mf.Seizure_Duration(Seizure_Count(1,i_ch_out),i_ch_out) = Seizure_Duration(Seizure_Count(1,i_ch_out),i_ch_out);
         mf.duration_amp_freq(Seizure_Count(1,i_ch_out),(i_ch_out-1)*3+1) = Seizure_Duration(Seizure_Count(1,i_ch_out),i_ch_out);
         
         %Also save these variables to the individual files (CKM)
         if DoSaveData
            %TODO: combine into single save if too much time?
            %Check if last spike occurs in next file. In this case Seizure_Count_Current is off, because we just reset it to 0s.
            mf2.last_spike_time_Current(Seizure_Count_Current(1,i_ch_out),i_ch_out) = last_spike_time(1,i_ch_out);
            mf2.Seizure_Duration_Current(Seizure_Count_Current(1,i_ch_out),i_ch_out) = Seizure_Duration(Seizure_Count(1,i_ch_out),i_ch_out);
         end
         
        
%         % make stim / no stim histogram
%         amp_vec = str2num(get(handles.ET_ampLo,'String'));
%         stim_on_logical = duration_amp_freq(:,(i_ch_out-1)*3+2) == amp_vec(i_ch_out);
%         stim_on_rows = duration_amp_freq(stim_on_logical,(i_ch_out-1)*3+(1:3));
%         stim_off_rows = duration_amp_freq(not(stim_on_logical),(i_ch_out-1)*3+(1:3));
%         
%         stim_off_rows(stim_off_rows(:,1)==0,:)=[]; % remove rows with zero seizure duration
%         
%         stim_on_durs = stim_on_rows(:,1);
%         stim_off_durs = stim_off_rows(:,1);
%         figure(i_ch_out)
%         subplot(1,2,1)
%         histogram(stim_on_durs,.5:1:90)
%         ylim([0 100])
%         
%         title(['i_ch_out = ' num2str(i_ch_out) ' stim'])
%         subplot(1,2,2)
%         histogram(stim_off_durs,.5:1:90)
%         title('no stim')
%         ylim([0 100])
        
        
        %% Bayes Opt for next stimulation parameters
        InitialObjective = log(Seizure_Duration(1:Seizure_Count(1,i_ch_out),i_ch_out));
        freq = stim_freq(1:Seizure_Count(1,i_ch_out),i_ch_out);
        amp = stim_amp(1:Seizure_Count(1,i_ch_out),i_ch_out);
        wid = stim_width(1:Seizure_Count(1,i_ch_out),i_ch_out);
        ampRatio = stim_ampRatio(1:Seizure_Count(1,i_ch_out),i_ch_out);
        chgRatio = stim_chgRatio(1:Seizure_Count(1,i_ch_out),i_ch_out);
        
        
        %CKM: If we are continuing an old recording, then we should start with all those parameters
        %EFF TODO: maybe not recalculating on every step? then again, does not show up in profiler
        if prior_Seizure_Count(i_ch_out)>0
            InitialObjective = [log(prior_Seizure_Duration(1:prior_Seizure_Count(1,i_ch_out),i_ch_out)); InitialObjective];
            freq = [prior_stim_freq(1:prior_Seizure_Count(1,i_ch_out),i_ch_out); freq];
            amp = [prior_stim_amp(1:prior_Seizure_Count(1,i_ch_out),i_ch_out); amp];
            wid = [prior_stim_width(1:prior_Seizure_Count(1,i_ch_out),i_ch_out); wid];
            ampRatio = [prior_stim_ampRatio(1:prior_Seizure_Count(1,i_ch_out),i_ch_out); ampRatio];
            chgRatio = [prior_stim_chgRatio(1:prior_Seizure_Count(1,i_ch_out),i_ch_out); chgRatio];
        end
        
%         freq_val_vec = linspace(freqLo_vec(i_ch), freqHi_vec(i_ch), freqBins)';
        freq_val_vec = logspace(log10(freqLo_vec(i_ch)), log10(freqHi_vec(i_ch)), freqBins)';
        amp_val_vec = linspace(ampLo_vec(i_ch), ampHi_vec(i_ch), ampBins)';
        wid_val_vec = linspace(widLo_vec(i_ch), widHi_vec(i_ch), widBins)';
        ampRatio_val_vec = linspace(ampRatioLo_vec(i_ch), ampRatioHi_vec(i_ch), ampRatioBins)';
        chgRatio_val_vec = linspace(chgRatioLo_vec(i_ch), chgRatioHi_vec(i_ch), chgRatioBins)';
        
        freq_int_vec = (1:freqBins)';
        amp_int_vec = (1:ampBins)';
        wid_int_vec = (1:widBins)';
        ampRatio_int_vec = (1:ampRatioBins)';
        chgRatio_int_vec = (1:chgRatioBins)';
        
        freq_binned_idx = knnsearch(freq_val_vec,freq,'K',1);
        amp_binned_idx = knnsearch(amp_val_vec,amp,'K',1);
        wid_binned_idx = knnsearch(wid_val_vec,wid,'K',1);
        ampRatio_binned_idx = knnsearch(ampRatio_val_vec,ampRatio,'K',1);
        chgRatio_binned_idx = knnsearch(chgRatio_val_vec,chgRatio,'K',1);
        
        freq_binned = freq_int_vec(freq_binned_idx);
        amp_binned = amp_int_vec(amp_binned_idx);
        wid_binned = wid_int_vec(wid_binned_idx);
        ampRatio_binned = ampRatio_int_vec(ampRatio_binned_idx);
        chgRatio_binned = chgRatio_int_vec(chgRatio_binned_idx);
        
        InitialX = table(freq_binned, amp_binned, wid_binned, ampRatio_binned, chgRatio_binned);
        
        freq_range = [1 freqBins];
        amp_range = [1 ampBins];
        wid_range = [1 widBins];
        ampRatio_range = [1 ampRatioBins];
        chgRatio_range = [1 chgRatioBins];
        
        opt_freq = optimizableVariable('frequency',freq_range,'Type','integer','Optimize',optimize_vec(1)); % Hz
        opt_amp = optimizableVariable('amplitude',amp_range,'Type','integer','Optimize',optimize_vec(2)); % Volts
        opt_wid = optimizableVariable('width',wid_range,'Type','integer','Optimize',optimize_vec(3)); % ms
        opt_ampRatio = optimizableVariable('amplitudeRatio',ampRatio_range,'Type','integer','Optimize',optimize_vec(4));
        opt_chgRatio = optimizableVariable('chgRatio',chgRatio_range,'Type','integer','Optimize',optimize_vec(5));
        
        %% remove non-optimized dims from xtrace
        var_names = InitialX.Properties.VariableNames;
        for i_dim = 1:5
            if optimize_vec(i_dim) == false
                InitialX = removevars(InitialX,var_names{1,i_dim});
            end
        end
        
        % instead: create header then send the data, if necessary, use
        % parfeval on a single worker to send the data
        % header will include bytes: ch, sz#, variable#, ro, co, data_bytes)
        % header is followed by data
        % this could be done with a single parfeval which sends the data then listens for a response
        % the parfeval function requires a header, the data
        
        bayes_or_AB_str = get(handles.ET_bayes_or_AB,'String');
        
        if strcmpi(bayes_or_AB_str,'Bayes') == 1
            sz_num = Seizure_Count(1,i_ch_out);
            ready_for_next_sz(i_ch_out) = false;
            disp(['ready_for_next_sz = ' num2str(ready_for_next_sz)]);
            %ready_for_next_sz = ready_for_next_sz %EFF21 - why this line?
            parfeval(@BO_wrapper, 0, opt_freq, opt_amp, opt_wid, opt_ampRatio, opt_chgRatio, InitialX, InitialObjective, sz_num, q{1,i_ch_out});
%             BO_wrapper(opt_freq, opt_amp, opt_wid, opt_ampRatio, opt_chgRatio, InitialX, InitialObjective, sz_num, q{1,i_ch_out});

        elseif strcmpi(bayes_or_AB_str,'Hier') == 1
            
            global hier_model_settings P_t_g theta_t_g off_set_est P_orig theta_orig
            
            sz_num = Seizure_Count(1,i_ch_out);
            ready_for_next_sz(i_ch_out) = false;
            disp(['ready_for_next_sz = ' num2str(ready_for_next_sz)]);
            %ready_for_next_sz = ready_for_next_sz
            
            legacy.opt_freq = opt_freq; % legacy variables used by BO_wrapper get passed through and later saved so that the old optimizer can be run in post as well if needed for comparison
            legacy.opt_amp = opt_amp;
            legacy.opt_wid = opt_wid;
            legacy.opt_ampRatio = opt_ampRatio;
            legacy.opt_chgRatio = opt_chgRatio;
            legacy.InitialX = InitialX;
            legacy.InitialObjective = InitialObjective;
            
            x = [freq_binned(end), amp_binned(end), wid_binned(end)];
            
            sz_dur = InitialObjective(end); % good

            parfeval(@Hier_wrapper, 0, x, P_t_g(:,:,i_ch_out), theta_t_g(:,i_ch_out), sz_dur, off_set_est(1,i_ch_out), P_orig(:,:,i_ch_out), theta_orig(:,i_ch_out), hier_model_settings, legacy, sz_num, q{1,i_ch_out});
%             h_dat = Hier_UNwrapper( x, P_t_g(:,:,i_ch_out), theta_t_g(:,i_ch_out), sz_dur, off_set_est(1,i_ch_out), P_orig(:,:,i_ch_out), theta_orig(:,i_ch_out), hier_model_settings, legacy, sz_num, q{1,i_ch_out});
%             disp('h_dat returned!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
%             send(q{1,i_ch_out},h_dat);

%%%%%%%%%%%%%%%%%   Direct: stim via ParamsA on every detection (no no-stim)   %%%%%%%%%%%%%%%%%%%%%%
            
        elseif strcmpi(bayes_or_AB_str,'Direct') == 1
            
%             InitialObjective = [log(prior_Seizure_Duration(1:prior_Seizure_Count(1,i_ch_out),i_ch_out)); InitialObjective];
%             freq = [prior_stim_freq(1:prior_Seizure_Count(1,i_ch_out),i_ch_out); freq];
%             amp = [prior_stim_amp(1:prior_Seizure_Count(1,i_ch_out),i_ch_out); amp];
%             wid = [prior_stim_width(1:prior_Seizure_Count(1,i_ch_out),i_ch_out); wid];
%             ampRatio = [prior_stim_ampRatio(1:prior_Seizure_Count(1,i_ch_out),i_ch_out); ampRatio];
%             chgRatio 
            
            AB_mat_name = [handles.save_folder filesep() 'res_ch' num2str(i_ch_out) '_sz' num2str(Seizure_Count(i_ch_out)) '_nr' num2str(n_read) '.mat'];
            save(AB_mat_name,'InitialX','InitialObjective','freq','amp','wid','ampRatio','chgRatio','Seizure_Count','i_ch_out','-nocompression')  
            
            coin = 1;
            
            if coin == 1 % best point
                if i_ch_out == 1
                    NextPoint_val_full = str2num(get(handles.ET_non_optimized_val_vec_ch1,'String'));
                elseif i_ch_out == 2
                    NextPoint_val_full = str2num(get(handles.ET_non_optimized_val_vec_ch2,'String'));
                elseif i_ch_out == 3
                    NextPoint_val_full = str2num(get(handles.ET_non_optimized_val_vec_ch3,'String'));
                elseif i_ch_out == 4
                    NextPoint_val_full = str2num(get(handles.ET_non_optimized_val_vec_ch4,'String'));
                end
            else
                error('coin must be 1')
            end
            
            next_freq(1,i_ch_out) = NextPoint_val_full(1,1);
            next_amp(1,i_ch_out) = NextPoint_val_full(1,2);
            next_width(1,i_ch_out) = NextPoint_val_full(1,3);
            next_ampRatio(1,i_ch_out) = NextPoint_val_full(1,4); % this was wid_val_vec
            next_chgRatio(1,i_ch_out) = NextPoint_val_full(1,5); % this was wid_val_vec

            set(handles.T_NextFreq,'String',num2str(next_freq));
            set(handles.T_NextAmp,'String',num2str(next_amp));
            set(handles.T_NextWidth,'String',num2str(next_width));
            set(handles.T_NextAmpRatio,'String',num2str(next_ampRatio));
            set(handles.T_NextChgRatio,'String',num2str(next_chgRatio));
            
            %%%%%%%%%%%%

%%%%%%%%%%%%%%%%%   AvsB   %%%%%%%%%%%%%%%%%%%%%%
            
        elseif strcmpi(bayes_or_AB_str,'AvsB') == 1
            
%                         InitialObjective = [log(prior_Seizure_Duration(1:prior_Seizure_Count(1,i_ch_out),i_ch_out)); InitialObjective];
%             freq = [prior_stim_freq(1:prior_Seizure_Count(1,i_ch_out),i_ch_out); freq];
%             amp = [prior_stim_amp(1:prior_Seizure_Count(1,i_ch_out),i_ch_out); amp];
%             wid = [prior_stim_width(1:prior_Seizure_Count(1,i_ch_out),i_ch_out); wid];
%             ampRatio = [prior_stim_ampRatio(1:prior_Seizure_Count(1,i_ch_out),i_ch_out); ampRatio];
%             chgRatio 
            
            AB_mat_name = [handles.save_folder filesep() 'res_ch' num2str(i_ch_out) '_sz' num2str(Seizure_Count(i_ch_out)) '_nr' num2str(n_read) '.mat'];
            save(AB_mat_name,'InitialX','InitialObjective','freq','amp','wid','ampRatio','chgRatio','Seizure_Count','i_ch_out','-nocompression')  
            
            coin = randi(2);
            
            if coin == 1 % best point
                if i_ch_out == 1
                    NextPoint_val_full = str2num(get(handles.ET_non_optimized_val_vec_ch1,'String'));
                elseif i_ch_out == 2
                    NextPoint_val_full = str2num(get(handles.ET_non_optimized_val_vec_ch2,'String'));
                elseif i_ch_out == 3
                    NextPoint_val_full = str2num(get(handles.ET_non_optimized_val_vec_ch3,'String'));
                elseif i_ch_out == 4
                    NextPoint_val_full = str2num(get(handles.ET_non_optimized_val_vec_ch4,'String'));
                end
            elseif coin == 2 %  no stim control
                NextPoint_val_full = [100 0 .1 3 1]; % set the amplitude to zero
            else
                error('coin must be 1 or 2')
            end
            
            next_freq(1,i_ch_out) = NextPoint_val_full(1,1);
            next_amp(1,i_ch_out) = NextPoint_val_full(1,2);
            next_width(1,i_ch_out) = NextPoint_val_full(1,3);
            next_ampRatio(1,i_ch_out) = NextPoint_val_full(1,4); % this was wid_val_vec
            next_chgRatio(1,i_ch_out) = NextPoint_val_full(1,5); % this was wid_val_vec

            set(handles.T_NextFreq,'String',num2str(next_freq));
            set(handles.T_NextAmp,'String',num2str(next_amp));
            set(handles.T_NextWidth,'String',num2str(next_width));
            set(handles.T_NextAmpRatio,'String',num2str(next_ampRatio));
            set(handles.T_NextChgRatio,'String',num2str(next_chgRatio));
            
            %%%%%%%%%%%%
            
            

%%%%%%%%%%%%%%%%%   ABC   %%%%%%%%%%%%%%%%%%%%%%
            
        elseif strcmpi(bayes_or_AB_str,'ABC') == 1
            
%                         InitialObjective = [log(prior_Seizure_Duration(1:prior_Seizure_Count(1,i_ch_out),i_ch_out)); InitialObjective];
%             freq = [prior_stim_freq(1:prior_Seizure_Count(1,i_ch_out),i_ch_out); freq];
%             amp = [prior_stim_amp(1:prior_Seizure_Count(1,i_ch_out),i_ch_out); amp];
%             wid = [prior_stim_width(1:prior_Seizure_Count(1,i_ch_out),i_ch_out); wid];
%             ampRatio = [prior_stim_ampRatio(1:prior_Seizure_Count(1,i_ch_out),i_ch_out); ampRatio];
%             chgRatio 
            
            AB_mat_name = [handles.save_folder filesep() 'res_ch' num2str(i_ch_out) '_sz' num2str(Seizure_Count(i_ch_out)) '_nr' num2str(n_read) '.mat'];
            save(AB_mat_name,'InitialX','InitialObjective','freq','amp','wid','ampRatio','chgRatio','Seizure_Count','i_ch_out','-nocompression')  
            
            coin = randi(3);
            
            if coin == 1 % best point
                if i_ch_out == 1
                    NextPoint_val_full = str2num(get(handles.ET_non_optimized_val_vec_ch1,'String'));
                elseif i_ch_out == 2
                    NextPoint_val_full = str2num(get(handles.ET_non_optimized_val_vec_ch2,'String'));
                elseif i_ch_out == 3
                    NextPoint_val_full = str2num(get(handles.ET_non_optimized_val_vec_ch3,'String'));
                elseif i_ch_out == 4
                    NextPoint_val_full = str2num(get(handles.ET_non_optimized_val_vec_ch4,'String'));
                end
            elseif coin ==2 % worst point
                if i_ch_out == 1
                    NextPoint_val_full = str2num(get(handles.ET_non_optimized_val_vec_ch1_worst,'String'));
                elseif i_ch_out == 2
                    NextPoint_val_full = str2num(get(handles.ET_non_optimized_val_vec_ch2_worst,'String'));
                elseif i_ch_out == 3
                    NextPoint_val_full = str2num(get(handles.ET_non_optimized_val_vec_ch3_worst,'String'));
                elseif i_ch_out == 4
                    NextPoint_val_full = str2num(get(handles.ET_non_optimized_val_vec_ch4_worst,'String'));
                end
            elseif coin == 3 %  no stim control
                NextPoint_val_full = [100 0 .1 3 1]; % set the amplitude to zero
            else
                error('coin must be 1, 2 or 3')
            end
            
            next_freq(1,i_ch_out) = NextPoint_val_full(1,1);
            next_amp(1,i_ch_out) = NextPoint_val_full(1,2);
            next_width(1,i_ch_out) = NextPoint_val_full(1,3);
            next_ampRatio(1,i_ch_out) = NextPoint_val_full(1,4); % this was wid_val_vec
            next_chgRatio(1,i_ch_out) = NextPoint_val_full(1,5); % this was wid_val_vec

            set(handles.T_NextFreq,'String',num2str(next_freq));
            set(handles.T_NextAmp,'String',num2str(next_amp));
            set(handles.T_NextWidth,'String',num2str(next_width));
            set(handles.T_NextAmpRatio,'String',num2str(next_ampRatio));
            set(handles.T_NextChgRatio,'String',num2str(next_chgRatio));
            
            %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%     ABCD     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
                  
        elseif strcmpi(bayes_or_AB_str,'ABCD') == 1
            
%                         InitialObjective = [log(prior_Seizure_Duration(1:prior_Seizure_Count(1,i_ch_out),i_ch_out)); InitialObjective];
%             freq = [prior_stim_freq(1:prior_Seizure_Count(1,i_ch_out),i_ch_out); freq];
%             amp = [prior_stim_amp(1:prior_Seizure_Count(1,i_ch_out),i_ch_out); amp];
%             wid = [prior_stim_width(1:prior_Seizure_Count(1,i_ch_out),i_ch_out); wid];
%             ampRatio = [prior_stim_ampRatio(1:prior_Seizure_Count(1,i_ch_out),i_ch_out); ampRatio];
%             chgRatio 
            
            AB_mat_name = [handles.save_folder filesep() 'res_ch' num2str(i_ch_out) '_sz' num2str(Seizure_Count(i_ch_out)) '_nr' num2str(n_read) '.mat'];
            save(AB_mat_name,'InitialX','InitialObjective','freq','amp','wid','ampRatio','chgRatio','Seizure_Count','i_ch_out','-nocompression')  
            
            coin = randi(4);
            
            if coin == 1 % ParamA
                if i_ch_out == 1
                    NextPoint_val_full = str2num(get(handles.ET_non_optimized_val_vec_ch1,'String'));
                elseif i_ch_out == 2
                    NextPoint_val_full = str2num(get(handles.ET_non_optimized_val_vec_ch2,'String'));
                elseif i_ch_out == 3
                    NextPoint_val_full = str2num(get(handles.ET_non_optimized_val_vec_ch3,'String'));
                elseif i_ch_out == 4
                    NextPoint_val_full = str2num(get(handles.ET_non_optimized_val_vec_ch4,'String'));
                end
            elseif coin ==2 % ParamB
                if i_ch_out == 1
                    NextPoint_val_full = str2num(get(handles.ET_non_optimized_val_vec_ch1_worst,'String'));
                elseif i_ch_out == 2
                    NextPoint_val_full = str2num(get(handles.ET_non_optimized_val_vec_ch2_worst,'String'));
                elseif i_ch_out == 3
                    NextPoint_val_full = str2num(get(handles.ET_non_optimized_val_vec_ch3_worst,'String'));
                elseif i_ch_out == 4
                    NextPoint_val_full = str2num(get(handles.ET_non_optimized_val_vec_ch4_worst,'String'));
                end
           elseif coin ==3 % ParamC
                if i_ch_out == 1
                    NextPoint_val_full = str2num(get(handles.ET_ParamC_ch1,'String'));
                elseif i_ch_out == 2
                    NextPoint_val_full = str2num(get(handles.ET_ParamC_ch2,'String'));
                elseif i_ch_out == 3
                    NextPoint_val_full = str2num(get(handles.ET_ParamC_ch3,'String'));
                elseif i_ch_out == 4
                    NextPoint_val_full = str2num(get(handles.ET_ParamC_ch4,'String'));
                end
            elseif coin == 4 %  no stim control
                NextPoint_val_full = [100 0 .1 3 1]; % set the amplitude to zero
            else
                error('coin must be 1, 2, 3 or 4')
            end
            
            next_freq(1,i_ch_out) = NextPoint_val_full(1,1);
            next_amp(1,i_ch_out) = NextPoint_val_full(1,2);
            next_width(1,i_ch_out) = NextPoint_val_full(1,3);
            next_ampRatio(1,i_ch_out) = NextPoint_val_full(1,4); % this was wid_val_vec
            next_chgRatio(1,i_ch_out) = NextPoint_val_full(1,5); % this was wid_val_vec

            set(handles.T_NextFreq,'String',num2str(next_freq));
            set(handles.T_NextAmp,'String',num2str(next_amp));
            set(handles.T_NextWidth,'String',num2str(next_width));
            set(handles.T_NextAmpRatio,'String',num2str(next_ampRatio));
            set(handles.T_NextChgRatio,'String',num2str(next_chgRatio));
            
        else
            error('must be Bayes or ABC or ABCD')
        end

    end
end

disp(['time_out = ' num2str(time_out)])

for i_ch_out = 1:n_ch_out
    [bo_dat, gotMsg] = poll(q{1,i_ch_out}, .001);
    
    if strcmpi(bayes_or_AB_str,'Bayes') == 1
        if gotMsg
            
            res = bo_dat.res;

            %% save res
            res_mat_name = [handles.save_folder filesep() 'res_ch' num2str(i_ch_out) '_sz' num2str(bo_dat.sz_num) '_nr' num2str(n_read) '.mat'];
%             res_mat_name = [handles.save_folder filesep() 'res_ch' num2str(i_ch_out) '_nr' num2str(n_read) '.mat'];

            save(res_mat_name,'res','-v7.3','-nocompression')        %TODO; RAN ERROR WITHOUT '-V7.3'

            %% set the next params
            non_optimized_val_vec = str2num(get(handles.ET_non_optimized_val_vec,'String'));

            freq_val_vec = logspace(log10(freqLo_vec(i_ch)), log10(freqHi_vec(i_ch)), freqBins)';
            amp_val_vec = linspace(ampLo_vec(i_ch), ampHi_vec(i_ch), ampBins)';
            wid_val_vec = linspace(widLo_vec(i_ch), widHi_vec(i_ch), widBins)';
            ampRatio_val_vec = linspace(ampRatioLo_vec(i_ch), ampRatioHi_vec(i_ch), ampRatioBins)';
            chgRatio_val_vec = linspace(chgRatioLo_vec(i_ch), chgRatioHi_vec(i_ch), chgRatioBins)';

            NextPoint_dimsel = table2array(res.NextPoint);
            NextPoint_ind_full = ones(1,5); % defaults to index 1, but all should be overwritten
            NextPoint_ind_full(optimize_vec) = NextPoint_dimsel;

            NextPoint_val_full = zeros(1,5);
            NextPoint_val_full(1,1) = freq_val_vec(NextPoint_ind_full(1,1));
            NextPoint_val_full(1,2) = amp_val_vec(NextPoint_ind_full(1,2));
            NextPoint_val_full(1,3) = wid_val_vec(NextPoint_ind_full(1,3));
            NextPoint_val_full(1,4) = ampRatio_val_vec(NextPoint_ind_full(1,4));
            NextPoint_val_full(1,5) = chgRatio_val_vec(NextPoint_ind_full(1,5));

            NextPoint_val_full(not(optimize_vec)) = non_optimized_val_vec(not(optimize_vec));

            next_freq(1,i_ch_out) = NextPoint_val_full(1,1);
            next_amp(1,i_ch_out) = NextPoint_val_full(1,2);
            next_width(1,i_ch_out) = NextPoint_val_full(1,3);
            next_ampRatio(1,i_ch_out) = NextPoint_val_full(1,4); % this was wid_val_vec
            next_chgRatio(1,i_ch_out) = NextPoint_val_full(1,5); % this was wid_val_vec

            set(handles.T_NextFreq,'String',num2str(next_freq));
            set(handles.T_NextAmp,'String',num2str(next_amp));
            set(handles.T_NextWidth,'String',num2str(next_width));
            set(handles.T_NextAmpRatio,'String',num2str(next_ampRatio));
            set(handles.T_NextChgRatio,'String',num2str(next_chgRatio));
            
            ready_for_next_sz(i_ch_out) = true;
            disp(['ready_for_next_sz = ' num2str(ready_for_next_sz)]);
            %ready_for_next_sz = ready_for_next_sz
        end
    elseif strcmpi(bayes_or_AB_str,'Hier')
        if gotMsg
            
            hier_mat_name = [handles.save_folder filesep() 'hier_ch' num2str(i_ch_out) '_sz' num2str(bo_dat.sz_num) '_nr' num2str(n_read) '.mat'];

            save(hier_mat_name,'bo_dat','v7.3','-nocompression')
            
            theta_t_g(:,i_ch_out) = bo_dat.theta_tp1;
            P_t_g(:,:,i_ch_out) = bo_dat.P_tp1;
            off_set_est(1,i_ch_out) = bo_dat.off_set_est;
            
            
            NextPoint_dimsel = bo_dat.next_point;
            
            %% set the next params
            
            non_optimized_val_vec = str2num(get(handles.ET_non_optimized_val_vec,'String'));

            freq_val_vec = logspace(log10(freqLo_vec(i_ch)), log10(freqHi_vec(i_ch)), freqBins)';
            amp_val_vec = linspace(ampLo_vec(i_ch), ampHi_vec(i_ch), ampBins)';
            wid_val_vec = linspace(widLo_vec(i_ch), widHi_vec(i_ch), widBins)';
            ampRatio_val_vec = linspace(ampRatioLo_vec(i_ch), ampRatioHi_vec(i_ch), ampRatioBins)';
            chgRatio_val_vec = linspace(chgRatioLo_vec(i_ch), chgRatioHi_vec(i_ch), chgRatioBins)';

            NextPoint_ind_full = ones(1,5); % defaults to index 1, but all should be overwritten
            NextPoint_ind_full(optimize_vec) = NextPoint_dimsel;

            NextPoint_val_full = zeros(1,5);
            NextPoint_val_full(1,1) = freq_val_vec(NextPoint_ind_full(1,1));
            NextPoint_val_full(1,2) = amp_val_vec(NextPoint_ind_full(1,2));
            NextPoint_val_full(1,3) = wid_val_vec(NextPoint_ind_full(1,3));
            NextPoint_val_full(1,4) = ampRatio_val_vec(NextPoint_ind_full(1,4));
            NextPoint_val_full(1,5) = chgRatio_val_vec(NextPoint_ind_full(1,5));

            NextPoint_val_full(not(optimize_vec)) = non_optimized_val_vec(not(optimize_vec));

            next_freq(1,i_ch_out) = NextPoint_val_full(1,1);
            next_amp(1,i_ch_out) = NextPoint_val_full(1,2);
            next_width(1,i_ch_out) = NextPoint_val_full(1,3);
            next_ampRatio(1,i_ch_out) = NextPoint_val_full(1,4); % this was wid_val_vec
            next_chgRatio(1,i_ch_out) = NextPoint_val_full(1,5); % this was wid_val_vec

            set(handles.T_NextFreq,'String',num2str(next_freq));
            set(handles.T_NextAmp,'String',num2str(next_amp));
            set(handles.T_NextWidth,'String',num2str(next_width));
            set(handles.T_NextAmpRatio,'String',num2str(next_ampRatio));
            set(handles.T_NextChgRatio,'String',num2str(next_chgRatio));
            
            ready_for_next_sz(i_ch_out) = true;
            disp(['ready_for_next_sz = ' num2str(ready_for_next_sz)]);
            %ready_for_next_sz = ready_for_next_sz
            
        end
    end
end

% disp('Seizure_Duration = ')
% disp(num2str(Seizure_Duration))

% disp('duration_amp_freq = ')
% disp(num2str(duration_amp_freq))
end %end of ~DoExisting block
prev_digi_data = digi_data;
prev_digi_time = digi_time;


%% --- Called by Process_Plot_Save
function BO_wrapper(opt_freq, opt_amp, opt_wid, opt_ampRatio,  opt_chgRatio, InitialX, InitialObjective, sz_num, que)

bo_dat.res = bayesopt(@sin,[opt_freq, opt_amp, opt_wid, opt_ampRatio,  opt_chgRatio],'InitialX',InitialX,'InitialObjective',InitialObjective, 'MaxObjectiveEvaluations', 1, 'PlotFcn', [], 'AcquisitionFunctionName', 'expected-improvement-plus', 'ExplorationRatio', 0.5, 'GPActiveSetSize',250);
bo_dat.sz_num = sz_num;
pause(25) % wait 25 seconds to ensure the stim has ended before the next stim settings are set
send(que,bo_dat)

%% --- Executes using lh2 whenever more data is needed (buffer runs below 1.4s)
function Generate_Stim_Vec(src, event, handles)

global DoTroubleshoot
if DoTroubleshoot
   disp('Run: Generate_Stim_Vec (Called by listener lh2)...')
end

global stim_flag
global fs_rec
global out_chunk
global next_freq next_amp next_width next_ampRatio next_chgRatio
global Seizure_Count
global stim_remaining

n_ch_out = length(stim_flag);
dt = 1/fs_rec;
n_t = fix(fs_rec*out_chunk);
data=zeros(n_t, n_ch_out*2);

estim_or_opto = get(handles.ET_estim_or_opto,'String');

if strcmpi(estim_or_opto,'estim')
    train_dur_vec = str2num(get(handles.ET_TrainDuration,'String'));
    stim_remaining(stim_flag) = train_dur_vec(stim_flag);
elseif strcmpi(estim_or_opto,'opto')
    train_dur_vec = str2num(get(handles.T_NextAmpRatio,'String'));
    stim_remaining(stim_flag) = train_dur_vec(stim_flag);
else
    error('must be estim or opto')
end

% ratio = [1, 4]; % pulse width ratio of pulse1 to pulse2

% amp_ratio_vec = str2num(get(handles.ET_ampRatioLo,'String')); % this was wrong!
% reversibility_ratio__vec = str2num(get(handles.ET_chgRatioHi,'String'));
    
if any(stim_remaining>0) 
    for i_ch_out = 1:n_ch_out
        if stim_remaining(i_ch_out)>0
            
            if strcmpi(estim_or_opto,'estim')
                amp_ratio = [1  next_ampRatio(i_ch_out)];
                rev_ratio = [1  next_chgRatio(i_ch_out)];
                data(:,i_ch_out) = pulse_train(n_t,amp_ratio,rev_ratio,next_amp(i_ch_out),next_freq(i_ch_out), next_width(i_ch_out),stim_remaining(i_ch_out));
            elseif strcmpi(estim_or_opto,'opto')
                amp_ratio = [next_ampRatio(i_ch_out)];
                rev_ratio = [next_chgRatio(i_ch_out)];
                data(:,i_ch_out) = pulse_train_opto(n_t,amp_ratio,rev_ratio,next_amp(i_ch_out),next_freq(i_ch_out), next_width(i_ch_out), stim_remaining(i_ch_out));
            else
                error('must be esti or opto')
            end
            
            data(:,i_ch_out+n_ch_out) = 1; % train on/off digital
            
            t = (1:n_t)*dt;
%             plot_vec = 1:fix(fs_rec*.05);
%             plot(handles.A_NextStim,t(plot_vec),data((plot_vec ),i_ch_out)+(i_ch_out-1)*3,'b')
%             hold on
%             ylim(handles.A_NextStim,[-10 n_ch_out*10])
            
            stim_remaining(i_ch_out) = stim_remaining(i_ch_out)-out_chunk;
        end
    end
    hold off
    set(handles.T_StimRemaining,'String',num2str(stim_remaining));
end

% for i_ch_out = 1:n_ch_out
%     if stim_remaining(i_ch_out)<=0 % 
% 
%         coin = mod(Seizure_Count(1,i_ch_out),2);
% 
%         if coin == 1
%             freq_vec = str2num(get(handles.ET_freqLo,'String'));
%             next_freq(1,i_ch_out) = freq_vec(i_ch_out);
%             amp_vec = str2num(get(handles.ET_ampLo,'String'));
%             next_amp(1,i_ch_out) = amp_vec(i_ch_out);
%             width_vec = str2num(get(handles.ET_widthLo,'String'));
%             next_width(1,i_ch_out) = width_vec(i_ch_out);
%         else
%             next_freq(1,i_ch_out) = 10; % had to set it to something
%             next_amp(1,i_ch_out) = 0;
%             next_width(1,i_ch_out) = 0.123;
%         end
% 
%         set(handles.T_NextFreq,'String',num2str(next_freq));
%         set(handles.T_NextAmp,'String',num2str(next_amp));
%         set(handles.T_NextWidth,'String',num2str(next_width));
%     end
% end

queueOutputData(src,data);

%% --- Called by Generate_Stim_Vec
function [out] = pulse_train(n_t,ratio,rev_ratio,charge,frequency, width, stim_remaining)
% amplitude = uA --> nC
% frequency = Hz
% charge = nC
% width = ms

% rev_ratio(2) = 0 gives monophasic

global fs_rec out_chunk

width_seconds = width/1000;

width_samples = fix(width_seconds*fs_rec);

amplitude = charge*10^-9/width_seconds; % amps
amplitude = amplitude*10^3; % mA
amplitude = amplitude/.1; % 1V/.1mA 

% if biphasic
total_width = width_samples * sum(ratio(:));
pulse = zeros(total_width,1);
pulse(1:(width_samples*ratio(1))) = amplitude;
pulse((width_samples*ratio(1))+1:end) = -amplitude*(ratio(1)/ratio(2))*(rev_ratio(2)/rev_ratio(1));
% monophasic = rev_ratio = 0;

out = zeros(n_t,1);

period_samples = fix(fs_rec*1/frequency);

if stim_remaining<out_chunk % fraction of a second of stim
    out(fix(period_samples/2):period_samples:stim_remaining*fs_rec) = 1;
else
    out(fix(period_samples/2):period_samples:end) = 1;
end

% out(randi([-fix(.4*period_samples) fix(.4*period_samples)],n_t,1) + (fix(period_samples/2):period_samples:(n_t-fix(period_samples/2)))) = 1;\

out = filter(pulse,1,out);
out(end) = 0; % end on a zero for safety

%% --- Called by Generate_Stim_Vec
function [out] = pulse_train_opto(n_t,ratio,rev_ratio,charge,frequency, width, t_remaining)
% amplitude = uA --> nC
% frequency = Hz
% charge = nC
% width = ms

train_dur = ratio;
not_used = rev_ratio;
amplitude = charge;
frequency = frequency;
duty_cycle = width;

global fs_rec

pulse_width_samples = fix(duty_cycle*(fs_rec/frequency));
pulse = amplitude * ones(pulse_width_samples,1);

delta_pulses = zeros(fix(fs_rec*train_dur),1);

period_samples = fix(fs_rec*1/frequency);  % this doesn't deal with pulse trains that are fractional length very well

delta_pulses(1:period_samples:end) = 1;

% out(randi([-fix(.4*period_samples) fix(.4*period_samples)],n_t,1) + (fix(period_samples/2):period_samples:(n_t-fix(period_samples/2)))) = 1;\

full_out = filter(pulse,1,delta_pulses);

full_out(end) = 0;


if ceil(train_dur)>train_dur % if the train duration is fractional, pad with zeros up to a full second
    full_out = [full_out; zeros(fix(fs_rec*(ceil(train_dur)-train_dur)),1)];
end

if length(full_out) < fs_rec*train_dur
    error('output wrong length')
end

if any(full_out)>5
    error('out of range here!')
end
if any(full_out)<0
    error('out of range here!')
end

t_start = train_dur - t_remaining;
s_start = fix(t_start*fs_rec)+1;
t_end = t_start + 1;
s_end = fix(t_end*fs_rec);
out = full_out(s_start:s_end);

%% --- function is not currently used
function [z] = place_holder_fcn(x, y)
z = x+y; % not actually used.




%% HELPER FUNCTIONS FOR SPIKE DETECTION
%%=====================================================================

% --- SpikeFinder is called by ProcessPlotSave
function [posspikes_wide,negspikes_wide]=EffSpikeFinder(data,fs,Settings)
%OPTIMIZED BY REMOVING ALL VARS WE DONT USE IN CURRENT CODE
%Find the location and width of all spikes in the data that adhere to the spike-settings    
%
%Input: 'data'                Row vector containing data, sampled at rate 'fs'
%       'fs'                  Number of samples per second
%       'Settings'            Uses the components 'Settings.ht_perc' and
%                            'Settings.req_width'
%       'DoDisplay'           1 = show popup of spikes (different colors
%                                 for wide vs narrow spikes)
%                             0 = don't show popup
%
%Output: 'posspikes'          posspikes(i,:) [location, height, width]
%        'negspikes'          where location references 'data'.
%        'posspikes_narrow'   posspikes contains information on all
%        'negspikes_narrow'   positive spikes. Similarly negspikes
%        'posspikes_wide'     contains info on all negative spikes.
%        'negspikes_wide'     _narrow and _wide specify only narrow and
%                             wide spike information
%
%Calls:         'peekseak', 'WidthFinder'
%
%Called by:     'MenuAnalyze_Callback' in 'findallofthem.m'
%

%Find all spikes, using peekseak
box_pos = peakseek(data,  Settings.dist_min,Settings.SpikesUpper, Settings.SpikesPosMax, Settings.dist_max);        %positive peak loc & data
box_neg = peakseek(-data, Settings.dist_min,-1*Settings.SpikesLower, -1*Settings.SpikesNegMax, Settings.dist_max);     %negative peak loc & data

%Find the width and height information of all spikes. Also, separate the
%spikes into a list of _wide spikes and a list of _narrow spikes.
posspikes_wide=EffWidthFinder(box_pos,0,data,fs,Settings);      
negspikes_wide=EffWidthFinder(box_neg,1,data,fs,Settings);


function spikes_wide = EffWidthFinder(evnt, sign, data, fs, Settings)
%EFF21 - adapted to reduce parts we don't need
%  Given a vector of spikes (location and raw data), calculate the width of the spikes, measured at a certain user specified height.
%  Positive and negative spikes are considered separately and specified with the variable 'sign'
%
%Input:         'evnt'      Row vector of peak locations, stored as indices of 'data'
%               'sign'      indicates whether positive or negative spikes are considered
%                           0  = positive peaks
%                           1  = negative peaks
%               'data'      Row vector containing data, sampled at rate 'fs'
%               'fs'        Number of samples per second
%               'Settings'  Relevant settings components:
%                           'Settings.ht_perc'  : height of a spike at
%                           which to determine its width (fraction)
%                           'Settings.req_width': width cutoff between
%                           narrow and wide spikes (s)
%
%Output:        'spikelist'     List of ALL peaks (both narrow and wide)
%                               spikelist(i,:)=[location, height, width]
%                               where location references the index in 'data'
%               'spikes_wide'   List of only the 'wide' peak info
%               'spikes_narrow' List of only the 'narrow' peak info
%               'spikes_too_wide' List of only the 'too wide' peak info
%
%Calls:         Flip
%
%Called by:     SpikeFinder

%Load in relevant Settings
ht_perc=Settings.ht_perc;       %Width of a spike is determined at specified height of spikes (fraction)
req_width=Settings.req_width;   %Spike is considered wide if width exceeds this threshold (in s)

%Initialize structures
spikes_wide=[];

%Loop through all peaks found. Find the width of each peak. Sometimes
%the data will change direction before height threshold is reached in
%which case interpolation might be used to calculate the peak height.
for i = 1:length(evnt)
    pk = evnt(i);                                        %location of peak
    pk_ht = data(pk);                                    %height of peak
    pk_ht_perc = data(pk) * ht_perc;                     %height of pk thresh
    left = pk;  %will be used to find width
    right = pk; %will be used to find witth

    %First attempt to just walk down the graph to the height threshold
    %This will only work if the graph does not go up again before
    %reaching it within the specified width
    leftsuccess=1;
    while leftsuccess &&...
          Flip((pk_ht_perc<data(left)),sign)                
      if left==1
          leftsuccess=0;
      elseif pk-left<0.5*req_width*fs  %making sure it terminates
         left=left-1;
      else
          leftsuccess=0;
      end
    end
    if leftsuccess
        %We found left edge of spike
    else
        if left==1
            %Peak is too far to the beginning of file, we cannot locate
            %left of spike. Double the right portion for best estimate
        else
            %Too much meandering. Use extrapolation to find width
            left = floor(pk - abs((pk - left) * ...
            (pk_ht - pk_ht_perc) / (pk_ht - data(left))));   %based on ratio
            leftsuccess=1;
        end
    end

    rightsuccess=1;
    while rightsuccess &&...
          Flip((pk_ht_perc<data(right)),sign)

      if right==length(data)
          rightsuccess=0;
      elseif right-pk<0.5*req_width*fs
          right=right+1;
      else
          rightsuccess=0;
      end
    end

    if rightsuccess
       %Right is the correct location 
    else
        if right==length(data)
            %Peak is too far to the end of file, we cannot locate right
            %of spike. Double left portion for best estimate
        else
            %Too much meandering. Use extrapolation to find width
            right = ceil(pk + abs((right - pk) * ...
            (pk_ht - pk_ht_perc) / (pk_ht - data(right))));   %based on ratio
            rightsuccess=1;
        end
    end

    %Deal with edge cases: 
    if ~leftsuccess && rightsuccess
        left=pk-(right-pk);
    elseif leftsuccess && ~rightsuccess
        right=pk+(pk-left);
    end

    %left and right are both successful, so calculate width
    spk_width = (right - left) / fs;

    if (spk_width<=Settings.req_width_max) && (spk_width>=Settings.req_width)
        spikes_wide=[spikes_wide ; [pk pk_ht spk_width]];
    end
end

function new_value = Flip(value, sign)
%Outputs boolean value based on 'sign'
%
%Input:     'value'     Boolean value (0 or 1)
%           'sign'      0 = return same value
%                       1 = return opposite value
%
%Output:    'new_value' boolean (0 or 1)
%
%Called by: WidthFinder
%
new_value = value;    
if sign == 1
   new_value = ~value;
end

function [locs, pks]=peakseek(x,minpeakdist,minpeakh,maxpeakh,maxpeakdist)
% Alternative to the findpeaks function.  This thing runs much much faster.
% It really leaves findpeaks in the dust.  It also can handle ties between
% peaks.  Findpeaks just erases both in a tie.  Shame on findpeaks.
%
% x is a vector input (generally a timecourse)
% minpeakdist is the minimum desired distance between peaks (optional, defaults to 1)
% minpeakh is the minimum height of a peak (optional)
%
% (c) 2010
% Peter O'Connor
% peter<dot>ed<dot>oconnor .AT. gmail<dot>com
%
% ***maxpeakh and maxpeakdist functionality added by Bethany Stieve***

if size(x,3)==1, x=x'; end

% Find all maxima and ties
locs=find(x(2:end-1)>=x(1:end-2) & x(2:end-1)>=x(3:end))+1;

if nargin<2, minpeakdist=1; end % If no minpeakdist specified, default to 1.

if nargin>2 % If there's a minpeakheight
    locs(x(locs)<=minpeakh)=[];
end

if nargin>3 % If there's a max peakheight
    locs(x(locs)> maxpeakh)=[];
end

if nargin<4, maxpeakdist=1; end % If no maxpeakdist specified, default to 1.

if minpeakdist>1
    while 1

        del=diff(locs)<minpeakdist;

        if ~any(del), break; end

        pks=x(locs);

        [garb mins]=min([pks(del) ; pks([false del])]); %#ok<ASGLU>

        deln=find(del);

        deln=[deln(mins==1) deln(mins==2)+1];

        locs(deln)=[];

    end
end

if maxpeakdist>1
    while 1

        del=diff(locs)>maxpeakdist;

        if ~any(del), break; end

        pks=x(locs);

        [garb mins]=min([pks(del) ; pks([false del])]); %#ok<ASGLU>

        deln=find(del);

        deln=[deln(mins==1) deln(mins==2)+1];

        locs(deln)=[];

    end
end

if nargout>1
    pks=x(locs);
end





%% CKM: FUNCTIONS FOR HIERARCHY
% =========================================================================

% --- Hier_wrapper is called by ProcessPlotSave
function Hier_wrapper(x, P_t, theta_t, sz_dur, off_set_est, P_orig, theta_orig, hier_model_settings, legacy, sz_num, que)

H_tp1 = zeros(1,hier_model_settings.n_setting_combos*(hier_model_settings.n_prior_mice+2)); % preallocate H_tp1
i_mouse = hier_model_settings.n_prior_mice+1; % the index of the new mouse is always the next mouse after the prior mice
h_ind = (i_mouse-1)*hier_model_settings.n_setting_combos+sub2ind(hier_model_settings.bins_of_optimized_dims, x(1), x(2), x(3));
H_tp1( h_ind ) = 1;

r_tp1 = sz_dur-off_set_est;

sigma = hier_model_settings.sigma;

G_t = P_t*H_tp1'*inv(H_tp1*P_t*H_tp1'+sigma^2);

theta_tp1 = theta_t + G_t*(r_tp1-H_tp1*theta_t);

P_tp1 = P_t-P_t*H_tp1'*inv(H_tp1*P_t*H_tp1'+sigma^2)*H_tp1*P_t;

%% check if mean is close enough to off_set_est, if not recalculate the entire thing
theta_subject = theta_tp1((i_mouse-1)*hier_model_settings.n_setting_combos+(1:hier_model_settings.n_setting_combos));
% disp('off_set below - - - - - - -  - - - - --- - - - - - - - - - - - - --')
off_set = mean(theta_subject(:))
off_set_est = off_set_est
% if and( abs(off_set)>0.15, sz_num > 10)
if and( abs(off_set)>0.1, sz_num > 3)
% if mod(sz_num,3)==2
    
    disp('recalculating !@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!')
    
    n_recalcs = 4;
    
    for i_recalc = 1:n_recalcs
        off_set = mean(theta_subject(:))
        off_set_est = off_set_est+off_set
 
        Hier_dat.recalculated = true;
        InitialX = table2array(legacy.InitialX);
        InitialObjective = legacy.InitialObjective;
        n_sz = length(InitialObjective);
        
        P_t = P_orig;
        theta_t = theta_orig;
        
        for i_sz = 1:n_sz

            H_tp1 = zeros(1,hier_model_settings.n_setting_combos*(hier_model_settings.n_prior_mice+2)); % preallocate H_tp1
            x_re = InitialX(i_sz,:);
            h_ind = (i_mouse-1)*hier_model_settings.n_setting_combos+sub2ind(hier_model_settings.bins_of_optimized_dims, x_re(1), x_re(2), x_re(3));
            H_tp1( h_ind ) = 1;

            r_tp1 = InitialObjective(i_sz)-off_set_est;

            sigma = hier_model_settings.sigma;

            G_t = P_t*H_tp1'*inv(H_tp1*P_t*H_tp1'+sigma^2);

            theta_tp1 = theta_t + G_t*(r_tp1-H_tp1*theta_t);

            P_tp1 = P_t-P_t*H_tp1'*inv(H_tp1*P_t*H_tp1'+sigma^2)*H_tp1*P_t;
            
            P_t = P_tp1;
            theta_t = theta_tp1;
        end
        
        theta_subject = theta_tp1((i_mouse-1)*hier_model_settings.n_setting_combos+(1:hier_model_settings.n_setting_combos)); % update the mean estimate
    end
else
    Hier_dat.recalculated = false;
end

%%
Hier_dat.off_set_est = off_set_est; % should update off_set_est global

Hier_dat.theta_tp1 = theta_tp1;
Hier_dat.P_tp1 = P_tp1;

Hier_dat.legacy = legacy;
Hier_dat.sz_num = sz_num;

% select next point using acqusition function !!!!

% Hier_dat.next_point = [1 2 3]; % for now use just 1 2 3 for testing

theta3D = reshape(theta_tp1((i_mouse-1)*hier_model_settings.n_setting_combos+(1:hier_model_settings.n_setting_combos)),hier_model_settings.bins_of_optimized_dims);
sig_diag = P_tp1(logical(eye(size(P_tp1))));     
sig3D = reshape(sig_diag((i_mouse-1)*hier_model_settings.n_setting_combos+(1:hier_model_settings.n_setting_combos)),hier_model_settings.bins_of_optimized_dims);
beta = 4;
[LCB3D, minX] = LCB(theta3D, sig3D, beta);
Hier_dat.LCB3D = LCB3D;
Hier_dat.next_point = minX;

send(que,Hier_dat)

%% unwrapper
function Hier_dat = Hier_UNwrapper(x, P_t, theta_t, sz_dur, off_set_est, P_orig, theta_orig, hier_model_settings, legacy, sz_num, que)

H_tp1 = zeros(1,hier_model_settings.n_setting_combos*(hier_model_settings.n_prior_mice+2)); % preallocate H_tp1
i_mouse = hier_model_settings.n_prior_mice+1; % the index of the new mouse is always the next mouse after the prior mice
h_ind = (i_mouse-1)*hier_model_settings.n_setting_combos+sub2ind(hier_model_settings.bins_of_optimized_dims, x(1), x(2), x(3));
H_tp1( h_ind ) = 1;

r_tp1 = sz_dur-off_set_est;

sigma = hier_model_settings.sigma;

G_t = P_t*H_tp1'*inv(H_tp1*P_t*H_tp1'+sigma^2);

theta_tp1 = theta_t + G_t*(r_tp1-H_tp1*theta_t);

P_tp1 = P_t-P_t*H_tp1'*inv(H_tp1*P_t*H_tp1'+sigma^2)*H_tp1*P_t;

%% check if mean is close enough to off_set_est, if not recalculate the entire thing
theta_subject = theta_tp1((i_mouse-1)*hier_model_settings.n_setting_combos+(1:hier_model_settings.n_setting_combos));
% disp('off_set below - - - - - - -  - - - - --- - - - - - - - - - - - - --')
off_set = mean(theta_subject(:))
off_set_est = off_set_est
% if and( abs(off_set)>0.15, sz_num > 10)
if and( abs(off_set)>0.05, sz_num > 2)
% if mod(sz_num,3)==2
    
    disp('recalculating !@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!@#!')
    
    n_recalcs = 4;
    
    for i_recalc = 1:n_recalcs
        off_set = mean(theta_subject(:))
        off_set_est = off_set_est+off_set
 
        Hier_dat.recalculated = true;
        InitialX = table2array(legacy.InitialX);
        InitialObjective = legacy.InitialObjective;
        n_sz = length(InitialObjective);
        
        P_t = P_orig;
        theta_t = theta_orig;
        
        for i_sz = 1:n_sz

            H_tp1 = zeros(1,hier_model_settings.n_setting_combos*(hier_model_settings.n_prior_mice+2)); % preallocate H_tp1
            x_re = InitialX(i_sz,:);
            h_ind = (i_mouse-1)*hier_model_settings.n_setting_combos+sub2ind(hier_model_settings.bins_of_optimized_dims, x_re(1), x_re(2), x_re(3));
            H_tp1( h_ind ) = 1;

            r_tp1 = InitialObjective(i_sz)-off_set_est;

            sigma = hier_model_settings.sigma;

            G_t = P_t*H_tp1'*inv(H_tp1*P_t*H_tp1'+sigma^2);

            theta_tp1 = theta_t + G_t*(r_tp1-H_tp1*theta_t);

            P_tp1 = P_t-P_t*H_tp1'*inv(H_tp1*P_t*H_tp1'+sigma^2)*H_tp1*P_t;
            
            P_t = P_tp1;
            theta_t = theta_tp1;
        end
        
        theta_subject = theta_tp1((i_mouse-1)*hier_model_settings.n_setting_combos+(1:hier_model_settings.n_setting_combos)); % update the mean estimate
    end
else
    Hier_dat.recalculated = false;
end

%%
Hier_dat.off_set_est = off_set_est; % should update off_set_est global

Hier_dat.theta_tp1 = theta_tp1;
Hier_dat.P_tp1 = P_tp1;

Hier_dat.legacy = legacy;
Hier_dat.sz_num = sz_num;

% select next point using acqusition function !!!!

% Hier_dat.next_point = [1 2 3]; % for now use just 1 2 3 for testing

theta3D = reshape(theta_tp1((i_mouse-1)*hier_model_settings.n_setting_combos+(1:hier_model_settings.n_setting_combos)),hier_model_settings.bins_of_optimized_dims);
sig_diag = P_tp1(logical(eye(size(P_tp1))));     
sig3D = reshape(sig_diag((i_mouse-1)*hier_model_settings.n_setting_combos+(1:hier_model_settings.n_setting_combos)),hier_model_settings.bins_of_optimized_dims);
beta = 4;
[LCB3D, minX] = LCB(theta3D, sig3D, beta);
Hier_dat.LCB3D = LCB3D;
Hier_dat.next_point = minX;

% send(que,Hier_dat)

function [LCB3D, minX] = LCB(theta3D, sig3D, beta)
    LCB3D = theta3D-beta*sig3D;
    
    [minM idx] = min(LCB3D(:));
    [d1 d2 d3] = ind2sub(size(LCB3D),idx);
    minX = [d1 d2 d3];

    
    
    
%% GUI FUNCTIONS FOR EACH BUTTON AND FIELD - NOT INTERESTING - NOTHING HAPPENS IN THEM
% =========================================================================
function ET_fs_Callback(hObject, eventdata, handles)
function ET_fs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_ampLo_Callback(hObject, eventdata, handles)
function ET_ampLo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_widthLo_Callback(hObject, eventdata, handles)
function ET_widthLo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_MonoOrBiPhasic_Callback(hObject, eventdata, handles)
function ET_MonoOrBiPhasic_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_freqLo_Callback(hObject, eventdata, handles)
function ET_freqLo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_TrainDuration_Callback(hObject, eventdata, handles)
function ET_TrainDuration_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_SaveFolder_Callback(hObject, eventdata, handles)
function ET_SaveFolder_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_SaveName_Callback(hObject, eventdata, handles)
function ET_SaveName_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_n_cams_Callback(hObject, eventdata, handles)
function ET_n_cams_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_n_ch_out_Callback(hObject, eventdata, handles)
function ET_n_ch_out_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_n_ch_in_Callback(hObject, eventdata, handles)
function ET_n_ch_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_seizure_detection_channels_Callback(hObject, eventdata, handles)
function ET_seizure_detection_channels_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_G_fast_Callback(hObject, eventdata, handles)
function ET_G_fast_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_G_slow_Callback(hObject, eventdata, handles)
function ET_G_slow_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_fast_slow_ratio_thresh_Callback(hObject, eventdata, handles)
function ET_fast_slow_ratio_thresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_open_loop_Callback(hObject, eventdata, handles)
function ET_open_loop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_channel_spacing_Callback(hObject, eventdata, handles)
function ET_channel_spacing_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_spike_dist_min_Callback(hObject, eventdata, handles)
function ET_spike_dist_min_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_spike_dist_max_Callback(hObject, eventdata, handles)
function ET_spike_dist_max_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_threshold_pos_Callback(hObject, eventdata, handles)
function ET_threshold_pos_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_threshold_neg_Callback(hObject, eventdata, handles)
function ET_threshold_neg_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_width_at_percent_height_Callback(hObject, eventdata, handles)
function ET_width_at_percent_height_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_min_width_Callback(hObject, eventdata, handles)
function T_min_width_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_min_spikes_per_two_s_Callback(hObject, eventdata, handles)
function T_min_spikes_per_two_s_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_device_ID_Callback(hObject, eventdata, handles)
function ET_device_ID_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_time_plot_Callback(hObject, eventdata, handles)
function ET_time_plot_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_end_max_spikes_Callback(hObject, eventdata, handles)
function ET_end_max_spikes_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_per_n_seconds_Callback(hObject, eventdata, handles)
function ET_per_n_seconds_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_exclusion_time_Callback(hObject, eventdata, handles)
function ET_exclusion_time_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_stim_dev_Callback(hObject, eventdata, handles)
function ET_stim_dev_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_fs_stim_Callback(hObject, eventdata, handles)
function ET_fs_stim_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_ampRatioLo_Callback(hObject, eventdata, handles)
function ET_ampRatioLo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_chgRatioHi_Callback(hObject, eventdata, handles)
function ET_chgRatioHi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_ampHi_Callback(hObject, eventdata, handles)
function ET_ampHi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_widthHi_Callback(hObject, eventdata, handles)
function ET_widthHi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_freqHi_Callback(hObject, eventdata, handles)
function ET_freqHi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function PB_load_for_restart_Callback(hObject, eventdata, handles)
function PB_restart_Callback(hObject, eventdata, handles)
function ET_path_for_restart_Callback(hObject, eventdata, handles)
function ET_path_for_restart_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_ampBins_Callback(hObject, eventdata, handles)
function ET_ampBins_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_widthBins_Callback(hObject, eventdata, handles)
function ET_widthBins_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_freqBins_Callback(hObject, eventdata, handles)
function ET_freqBins_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_ampRatioBins_Callback(hObject, eventdata, handles)
function ET_ampRatioBins_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_ampRatioHi_Callback(hObject, eventdata, handles)
function ET_ampRatioHi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_chgRatioLo_Callback(hObject, eventdata, handles)
function ET_chgRatioLo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_chgRatioBins_Callback(hObject, eventdata, handles)
function ET_chgRatioBins_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_optimize_vec_Callback(hObject, eventdata, handles)
function ET_optimize_vec_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ET_optimize_vec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_non_optimized_val_vec_Callback(hObject, eventdata, handles)
function ET_non_optimized_val_vec_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_non_optimized_val_vec_ch1_Callback(hObject, eventdata, handles)
function ET_non_optimized_val_vec_ch1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_non_optimized_val_vec_ch2_Callback(hObject, eventdata, handles)
function ET_non_optimized_val_vec_ch2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_non_optimized_val_vec_ch3_Callback(hObject, eventdata, handles)
function ET_non_optimized_val_vec_ch3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_non_optimized_val_vec_ch4_Callback(hObject, eventdata, handles)
function ET_non_optimized_val_vec_ch4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_bayes_or_AB_Callback(hObject, eventdata, handles)
function ET_bayes_or_AB_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_non_optimized_val_vec_ch1_worst_Callback(hObject, eventdata, handles)
function ET_non_optimized_val_vec_ch1_worst_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_non_optimized_val_vec_ch2_worst_Callback(hObject, eventdata, handles)
function ET_non_optimized_val_vec_ch2_worst_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_non_optimized_val_vec_ch3_worst_Callback(hObject, eventdata, handles)
function ET_non_optimized_val_vec_ch3_worst_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_non_optimized_val_vec_ch4_worst_Callback(hObject, eventdata, handles)
function ET_non_optimized_val_vec_ch4_worst_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_load_path_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_min_width_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_min_spikes_per_two_s_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_non_optmize_val_vec_Callback(hObject, eventdata, handles)
function ET_non_optmize_val_vec_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_hier_model_mat_Callback(hObject, eventdata, handles)
function ET_hier_model_mat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_spike_dist_min_raw_Callback(hObject, eventdata, handles)
function ET_spike_dist_min_raw_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_spike_dist_max_raw_Callback(hObject, eventdata, handles)
function ET_spike_dist_max_raw_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_threshold_pos_raw_Callback(hObject, eventdata, handles)
function ET_threshold_pos_raw_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_threshold_neg_raw_Callback(hObject, eventdata, handles)
function ET_threshold_neg_raw_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_width_at_percent_height_raw_Callback(hObject, eventdata, handles)
function ET_width_at_percent_height_raw_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_min_width_raw_Callback(hObject, eventdata, handles)
function ET_min_width_raw_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_and_or_spikes_Callback(hObject, eventdata, handles)
function ET_and_or_spikes_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_and_window_Callback(hObject, eventdata, handles)
function ET_and_window_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_spike_logic_1_Callback(hObject, eventdata, handles)
function ET_spike_logic_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_spike_logic_2_Callback(hObject, eventdata, handles)
function ET_spike_logic_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_spike_logic_3_Callback(hObject, eventdata, handles)
function ET_spike_logic_3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_spike_logic_4_Callback(hObject, eventdata, handles)
function ET_spike_logic_4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_min_spikes_per_two_s_stop_Callback(hObject, eventdata, handles)
function ET_min_spikes_per_two_s_stop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_estim_or_opto_Callback(hObject, eventdata, handles)
function ET_estim_or_opto_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_ParamC_ch1_Callback(hObject, eventdata, handles)
function ET_ParamC_ch1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_ParamC_ch2_Callback(hObject, eventdata, handles)
function ET_ParamC_ch2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_ParamC_ch3_Callback(hObject, eventdata, handles)
function ET_ParamC_ch3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_ParamC_ch4_Callback(hObject, eventdata, handles)
function ET_ParamC_ch4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_max_width_raw_Callback(hObject, eventdata, handles)
function ET_max_width_raw_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_max_width_Callback(hObject, eventdata, handles)
function ET_max_width_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_maxAmp_raw_pos_Callback(hObject, eventdata, handles)
function ET_maxAmp_raw_pos_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_maxAmp_raw_neg_Callback(hObject, eventdata, handles)
function ET_maxAmp_raw_neg_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_maxAmp_pos_Callback(hObject, eventdata, handles)
function ET_maxAmp_pos_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_maxAmp_neg_Callback(hObject, eventdata, handles)
function ET_maxAmp_neg_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ET_load_path_Callback(hObject, eventdata, handles)
function ET_save_path_Callback(hObject, eventdata, handles)
function ET_save_path_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function InputRefresh_Callback(hObject, eventdata, handles)
function InputRefresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function InputFileName_Callback(hObject, eventdata, handles)
function InputFileName_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function InputStartAt_Callback(hObject, eventdata, handles)
function InputStartAt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% POPUP PANEL to see graph tests
%==========================================================================
% --- Executes on button press in ButtonClose.
function ButtonClose_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonClose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.TextFS,'String','Sampling Rate');
set(handles.TextFN,'String','File name loaded');
set(handles.TextChIn,'String','#Channels In');
set(handles.TextChOut,'String','#Channels Out');
set(handles.PanelTest,'Visible','off');
set(handles.ButtonStart,'Enable','off');
set(handles.ButtonStop,'Enable','on');
set(handles.ButtonClose,'Enable','off');
set(handles.btn_open,'Enable','on');
set(handles.PB_Go,'Enable','on');
set(handles.PB_Stop,'Enable','on');
set(handles.TextStatus,'String','');
global DoRestartOffline
DoRestartOffline=0;

% --- Executes on button press in ButtonStart.
function ButtonStart_Callback(hObject, eventdata, handles)
%This reruns the analysis offline, so no need to resave the data files,
%though we DO have to save the overarching summary data file.
global ReanalysisFolder ReanalysisFiles ReanalysisIndex
global DoRestartOffline DoPauseOffline DoStopOffline DoExisting
global DoSaveToFile DoTroubleShoot
global mf save_mat_path save_mat_path2 spike_count
global fast_int slow_int Seizure_On Seizure_Off Seizure_Duration
global fast_slow_ratio_trigger
global spikes_trigger spikes_trigger_stop
global last_spike_time Seizure_Start_First_Spike_Time
global time_out stim_flag
global DoLog DoResaveDataFiles DoFinishedFile
global q
global AlreadyUpdated
global TimeStamp_postfix
if DoLog
    logname=['Log_Bayes_' 'rerun' '__' datestr(now,'mm_dd_yyyy__HH_MM') '.txt'];
    diary(logname);
end

if DoTroubleShoot
    disp('Start offline analysis...');
end

%Keep track of whether we already updated (UpdateFiles) for this current file. Don't want to do it twice.
AlreadyUpdated=0; %Did we already call UpdateFiles for current file?
DoFinishedFile=0; %Are we ready to call UpdateFiles through Process_Plot_Save

%reset plot
hold(handles.A_MainPlot, 'off')
plot(handles.A_MainPlot,0,0);
hold(handles.A_MainPlot, 'on')

%update GUI
DoStopOffline=0;
DoExisting=1;
ResetOfflineGUI('off','on',handles);
set(handles.PB_Go,'Enable','off');
set(handles.PB_Stop,'Enable','off');

%Allow user to run current file only, or all files in the file list
DoCurrentOnly=get(handles.CheckThisOnly,'value');
if DoCurrentOnly
   FilesTodo=ReanalysisIndex; 
else
   FilesTodo=1:size(ReanalysisFiles,2);
end

%Extract info on #channels & initialize for the first file
CurrentFile=FilesTodo(1);
ReanalysisPath=[ReanalysisFolder ReanalysisFiles{CurrentFile}];
InputMatfile=matfile(ReanalysisPath);
SizeData=size(InputMatfile,'data');

%Check if time exists as separate variable or is included as first column of data
listOfVariables = who('-file', ReanalysisPath);
TimeIncluded=0+~ismember('time', listOfVariables); % returns true if it is first column

%Extract necessary input from file
SizeNrChOut=size(InputMatfile,'detect_data');
SizeNrChOut=SizeNrChOut(2)-TimeIncluded;            %nr output channels
SizeNrChIn=SizeData(2)-SizeNrChOut-TimeIncluded;    %nr input channels
fs_rec=InputMatfile.fs_rec;                         %sampling rate
n_ch_in=SizeNrChIn;
n_ch_out=SizeNrChOut;

%Patch it so it can work in offline mode but with saving data files
if DoExisting && DoResaveDataFiles
  for i_ch = 1:n_ch_out
    q{1,i_ch} = parallel.pool.PollableDataQueue;
  end
  global prior_Seizure_Count prior_Seizure_Duration prior_stim_freq prior_stim_amp prior_stim_width prior_stim_ampRatio prior_stim_chgRatio
    prior_Seizure_Count = zeros(1,4);
    prior_Seizure_Duration = [];
    prior_stim_freq = [];
    prior_stim_amp = [];
    prior_stim_width = [];
    prior_stim_ampRatio = [];
    prior_stim_chgRatio = [];
end


%Pull these from settings in GUI - assume user will make sure they match ???TODO build in check
seizure_detection_ch_vec = str2num(get(handles.ET_seizure_detection_channels,'String'));% detect seizures on ACH0 and ACH2, these are BNC-2090 numbers
ch_in_vec = str2num(get(handles.ET_n_ch_in,'String')); % hard coded 1 input channel

%Initialize variables that are expected
InitializeVarsInput(n_ch_in,n_ch_out);

%Define output files. Initialize them. ???CKM
DoSaveToFile=get(handles.CheckSave,'value');
if DoSaveToFile
    %Load in global variables needed to save to file

    %Define extra variables to be stored in the files, to be used for
    %identification purposes.
    save_folder_base = get(handles.InputFileName,'String');
    save_name = get(handles.ET_SaveName,'String');
    save_folder = save_folder_base;
    handles.save_folder = save_folder;

    LiveRecording=0;                                    %LiveRecording=0 means offline reanalysis, LiveRecording=1 means new recording
    ReferenceDir='..\';                                 %This means the _d_ analysis files are in the parent directory
    ID=save_name;                                       %identifying ID (???for now this is the save_name chosen
    TimeStampString=datestr(datetime,'mm-dd-yyyy_HH-MM-SS');  %current time stamps. Included in EVERY file saved.
    TimeStamp_postfix=[' (' TimeStampString ')'];             %Add this string to each data file file name
    TimeStamp=now;
    SettingsName=['settings' TimeStamp_postfix '.mat']; %name of the settings file to be used
    IncludedDataFiles='';
    
    %Create subdirectory for running the reanalysis
    try
        mkdir(save_folder);
    catch        
        warndlg('Cannot create reanalysis subfolder');
    end
    
    %save settings to file
    save_settings(save_folder,SettingsName, handles);
    
    %Extract crucial info ch_in and ch_out from first data fileso we can initia;ize at this
    %point
    
    %save the overarching data file
    %EFF21 - remove Seizure_On
    save_mat_path = [save_folder '\' save_name  TimeStamp_postfix '.mat'];
    save(save_mat_path,'fast_int','slow_int','Seizure_Off','stim_flag','Seizure_Duration','spike_count','fast_slow_ratio_trigger','spikes_trigger','spikes_trigger_stop','Seizure_Start_First_Spike_Time','last_spike_time','time_out',...
                   'LiveRecording','IncludedDataFiles','ReferenceDir','ID','SettingsName','TimeStamp','TimeStampString','-v7.3','-nocompression')
    
    %Optionally, set it up for each data file
    if DoResaveDataFiles
        global mf2 save_mat_path2_base ro_data co_data ro_ard co_ard ro_dd co_dd datatype  d_num
        d_num = 1;
        
        %Data variables get initialized to zeros and saved to file
        time=zeros(1,1);                                 %EFF21: time as separate variable
        data=zeros(1,n_ch_in+n_ch_out);                  %EFF21: next n_ch_in columns are the input channels, and the last n_ch_out columns are the output channels
        art_rem_data = zeros(1,n_ch_out,datatype);       %EFF21
        art_rem_data_raw = zeros(1,n_ch_out,datatype);   %EFF21
        detect_data = zeros(1,n_ch_out,datatype);    %EFF21
        detect_data_raw = zeros(1,n_ch_out,datatype);%EFF21
        ro_data=1; co_data=size(data,2);                        %Keep track of row we are in
        ro_dd=1; co_dd=size(detect_data,2);                     %Keep track of row we are in
        ro_ard=1; co_ard=size(art_rem_data,2);             %Keep track of row we are in

        %Define extra variables so that individual files can store some of the data
        %as a safety precaution when the large file gets corrupted.
        global last_spike_time_Current Seizure_Start_First_Spike_Time_Current Seizure_Duration_Current stim_Current
        last_spike_time_Current =  zeros(1,n_ch_out);
        Seizure_Start_First_Spike_Time_Current = zeros(1,n_ch_out);
        Seizure_Duration_Current = zeros(1,n_ch_out);
        stim_Current = zeros(1,5*n_ch_out); %n_ch_out*[freq,amp,width,ampratio,chgratio]

        save_mat_path2_base = [save_folder '\' save_name '_d_'] % c
        save_mat_path2 = [save_mat_path2_base  num2str(d_num) TimeStamp_postfix '.mat'];

        save(save_mat_path2,'time','data','art_rem_data','art_rem_data_raw','detect_data','detect_data_raw','fs_rec',...
                        'ID','TimeStamp','last_spike_time_Current','Seizure_Start_First_Spike_Time_Current','Seizure_Duration_Current','stim_Current','-v7.3','-nocompression')
        mf2 = matfile(save_mat_path2,'Writable',true);
        clear time data art_rem_data detect_data art_rem_data_raw 
        
        global next_freq next_amp next_width  next_ampRatio next_chgRatio
        next_freq = 10*ones(1,n_ch_out); % initialize arbitrarily to 10 Hz
        next_amp = 5*ones(1,n_ch_out); % initialize arbitrarily to 1 unit
        next_width = 0.150*ones(1,n_ch_out); % initialize arbitrarily to 1 unit
        next_ampRatio = 1*ones(1,n_ch_out); % initialize to a 4:1 amplitude ratio
        next_chgRatio = 1*ones(1,n_ch_out); % initialize to a 1:1 chg ratio

        set(handles.T_NextFreq,'String',num2str(next_freq));
        set(handles.T_NextAmp,'String',num2str(next_amp));
        set(handles.T_NextWidth,'String',num2str(next_width));
        next_freq(1,:) = 10;
        next_amp(1,:) = 5;
        next_width(1,:) = 0.1;
        next_ampRatio(1,:) = 1;
        next_chgRatio(1,:) = 1;

        set(handles.T_NextFreq,'String',num2str(next_freq));
        set(handles.T_NextAmp,'String',num2str(next_amp));
        set(handles.T_NextWidth,'String',num2str(next_width));
        set(handles.T_NextAmpRatio,'String',num2str(next_ampRatio));
        set(handles.T_NextChgRatio,'String',num2str(next_chgRatio));

    end
    
    %We initialized these vars to large size for saving to file, but now
    %these can be reduced for internal use, where we only need to save most
    %recent 10 rows. 11 is HARDCODED as 1+10
    clear spike_count
    
    %create matfile for access
    mf = matfile(save_mat_path,'Writable',true);
end

%Loop through all files selected.
for i_file=1:size(FilesTodo,2)
    DoFinishedFile=(i_file>1); %Flag that is used in Process_Plot_Save to show we are ready to run UpdateFiles and setup next file
    AlreadyUpdated=0; %Flag to prevent UpdateFiles from happening more than once - should be reset each new file
    if DoStopOffline
        break
    end
    %Get reference to current file
    CurrentFile=FilesTodo(i_file);
    ReanalysisPath=[ReanalysisFolder ReanalysisFiles{CurrentFile}];
    disp(['Running file ' num2str(i_file) '/' num2str(size(FilesTodo,2)) ': ' ReanalysisFiles{CurrentFile}]);
    locs=strfind(ReanalysisFolder,'\');
    set(handles.TextFN,'String',['Loaded: ' ReanalysisFolder(1:3) '...' ReanalysisFolder(locs(end-1):end) ReanalysisFiles{CurrentFile}]);
    set(handles.ListOfFiles,'value',CurrentFile);
    set(handles.TextStatus,'String',['Analyzing file ' num2str(i_file) '/' num2str(size(FilesTodo,2)) ': ' ReanalysisFiles{CurrentFile}]);
    guidata(hObject,handles);
    
    %Optionally save to file
    if DoSaveToFile
        mf.IncludedDataFiles=[mf.IncludedDataFiles ',' ReanalysisFiles{CurrentFile}];
    end
    
    %reference the matfile
    InputMatfile=matfile(ReanalysisPath);
    SizeData=size(InputMatfile,'data');       
        
    %Chunk data and loop through the code as if we're first seeing it
    chunksize=fs_rec;
    nrchunks=floor(SizeData(1)/chunksize);
    RefreshTime=str2num(get(handles.InputRefresh,'String')); %how often to refresh plot

    %Determine what chunk to start in (matching user choice in GUI)
    if TimeIncluded
        starttime=InputMatfile.data(2,1);
    else
        starttime=InputMatfile.time(2,1);
    end
    chunk_startat=max(1,floor((str2num(get(handles.InputStartAt,'String'))-starttime)/(chunksize/fs_rec)));

    %EFF21: load whole file, then extract from that. Takes a bit more memory space, but a lot less time accessing.
        LoadedData=InputMatfile.data(2:SizeData(1),:);

        if TimeIncluded
            LoadedTime=LoadedData(:,1);  %extract time from first column
        else
            LoadedTime=InputMatfile.time;             %if stored separately
        end

    for i_data=chunk_startat:nrchunks
        while DoPauseOffline
          pause(1);
        end
        disp(['Loading datachunk ' num2str(i_data) ' / ' num2str(nrchunks) '...']);
        if ~DoStopOffline
            %Determine start and stop index of what part of variable to load in
            i_start=(i_data-1)*chunksize+1;
            i_stop=i_start+chunksize-1;

            %EFF21 Extract part of variables data and time - took this out of loop to
            %make faster at cost of more memroy use
            %LoadedData=InputMatfile.data(i_start:i_stop,TimeIncluded+1:SizeData(2)); %extract data
            %if TimeIncluded
            %    LoadedTime=InputMatfile.data(i_start:i_stop,TimeIncluded);  %extract time from first column
            %else
            %    LoadedTime=InputMatfile.time(i_start:i_stop,:);             %if stored separately
            %end

            %Do everything we would normally do with this data if it came from daq
            event.Data=LoadedData(i_start:i_stop,TimeIncluded+1:end); %extract data);
            event.TimeStamps=LoadedTime(i_start:i_stop,:);
            Process_Plot_Save(DoExisting,event, handles.A_MainPlot,ch_in_vec,seizure_detection_ch_vec,n_ch_in,n_ch_out,'', handles);
            DoFinishedFile=0; %Flag should be set to 0 as soon as we are done getting ready for next file (used in process_plot_save
            if mod(i_data,RefreshTime)==0 %only update graph every 10 chunks
                drawnow;
            end
        else
            break;
        end
        DoRestartOffline=1; %If we press this button again, it is a restart
    end
    
end %end loop through all files selected

%If process was stopped and there is still a file to be updated, do it now.
if ~AlreadyUpdated
    UpdateFiles(n_ch_out,(~DoExisting || DoResaveDataFiles));
end
    
%All done, so reset GUI
ResetOfflineGUI('on','off',handles)
disp('Finished offline analysis to end of file');
set(handles.TextStatus,'String','Done');

% --- Executes on button press in ButtonStop.
function ButtonStop_Callback(hObject, eventdata, handles)
%Make sure offline analysis stops after current block of data
global DoStopOffline DoPauseOffline DoSaveToFile AlreadyUpdated DoResaveDataFiles
global Seizure_Count
global mf mf2 Sz_On_Index Sz_Off_Index Sz_Index

DoPauseOffline=0;   %reset pause to 0 in case we were in pause while stopping
DoStopOffline=1;    %flag that we want to stop
ResetOfflineGUI('on','off',handles);

if DoSaveToFile && ~AlreadyUpdated
    ch_out_vec = str2num(get(handles.ET_n_ch_out,'String'));
    n_ch_out = length(ch_out_vec);
    UpdateFiles(n_ch_out,DoResaveDataFiles)
end

% --- Executes on button press in ButtonPause.
function ButtonPause_Callback(hObject, eventdata, handles)
%Pause running the offline analyis. Press again to resume
% hObject    handle to ButtonPause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DoPauseOffline
if DoPauseOffline
    set(handles.ButtonPause,'String','Pause');
    DoPauseOffline=0;
    set(handles.TextStatus,'String','Offline analysis resumed');
else
    set(handles.ButtonPause,'String','Resume');
    DoPauseOffline=1;
    set(handles.TextStatus,'String','Offline analysis paused');
end
   
% --- Executes on button press in ButtonSaveData.
function ButtonSaveData_Callback(hObject, eventdata, handles)
%Will save the data so far in file specified - not yet working
% hObject    handle to ButtonSaveData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Save the offline detection analysis to file, to be opened just like a live
%run. This button will save current run from memory, all at once




% --- Executes on button press in ButtonPreview: previews entire file
% loaded
function ButtonPreview_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonPreview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%reset plot
hold(handles.A_MainPlot, 'off')
plot(handles.A_MainPlot,0,0);
hold(handles.A_MainPlot, 'on')

global ReanalysisFolder ReanalysisFiles ReanalysisIndex

set(handles.TextStatus,'String',['Creating preview of file: ' ReanalysisFiles{ReanalysisIndex}]);
drawnow;

ReanalysisPath=[ReanalysisFolder ReanalysisFiles{ReanalysisIndex}];
%Create matfile for the file selected
listOfVariables = who('-file', ReanalysisPath);
TimeIncluded=0+~ismember('time', listOfVariables); % returns true if time is stored in first column of data

InputMatfile=matfile(ReanalysisPath);
SizeData=size(InputMatfile,'data');
SizeNrChOut=size(InputMatfile,'detect_data');
SizeNrChOut=SizeNrChOut(2)-TimeIncluded;            %nr output channels
SizeNrChIn=SizeData(2)-SizeNrChOut-TimeIncluded;    %nr input channels
n_ch_in=SizeNrChIn;
n_ch_out=SizeNrChOut;

%Start with entire file view in plot
raw_deci = 20;                               %plot is downsampled by this amount
if TimeIncluded
    %Time is the first column from data files
    raw_d = InputMatfile.data(:,2:SizeData(2));  
    timestamps= InputMatfile.data(:,1);
else
    %Time is stored as a separate variable
    raw_d = InputMatfile.data;  
    timestamps= InputMatfile.time;
end

raw_d = raw_d - mean(raw_d,1);               %remove mean for each column
n_ch=SizeData(2)-TimeIncluded;
channel_scaling = str2num(get(handles.ET_channel_spacing,'String'));
channel_spacing = [1:n_ch-n_ch_out 2*(1:n_ch_out)];
plot(handles.A_MainPlot, timestamps(1:raw_deci:SizeData(1),:), raw_d(1:raw_deci:SizeData(1),:)*diag(1./channel_scaling)+repmat(channel_spacing, length(timestamps(1:raw_deci:SizeData(1))),1),'LineWidth',1)
xlim([timestamps(2) timestamps(SizeData(1))]);
set(handles.TextStatus,'String',['Current preview: ' ReanalysisFiles{ReanalysisIndex}]);

function ResetOfflineGUI(change1,change2,handles)
%reset the buttons once start or stop is pressed or it is done
set(handles.ButtonStart,'Enable',change1);
set(handles.ButtonStop,'Enable',change2);
set(handles.ButtonClose,'Enable',change1);
set(handles.ButtonPreview,'Enable',change1);
set(handles.CheckThisOnly,'Enable',change1');
set(handles.CheckSave,'Enable',change1);
if get(handles.CheckSave,'value')
    set(handles.InputFileName,'Enable','on');
else
    set(handles.InputFileName,'Enable','off');
end

% --- Executes on button press in CheckSave.: change whether output dir is
% editable
function CheckSave_Callback(hObject, eventdata, handles)
% hObject    handle to CheckSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckSave
if get(handles.CheckSave,'Value')
    set(handles.InputFileName,'Enable','on');
else
    set(handles.InputFileName,'Enable','off');
end

function ListOfFiles_Callback(hObject, eventdata, handles)
% hObject    handle to ListOfFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ListOfFiles contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ListOfFiles
global ReanalysisFolder ReanalysisFiles ReanalysisIndex
ReanalysisIndex=get(handles.ListOfFiles,'value');
locs=strfind(ReanalysisFolder,'\');
set(handles.TextFN,'String',['Loaded: ' ReanalysisFolder(1:3) '...' ReanalysisFolder(locs(end-1):end) ReanalysisFiles{ReanalysisIndex}]);
guidata(hObject,handles);

function ListOfFiles_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CheckThisOnly_Callback(hObject, eventdata, handles)
function CheckFluid_Callback(hObject, eventdata, handles)


% --- Executes on button press in SelectSettings.

function SelectSettings_Callback(hObject, eventdata, handles)
% hObject    handle to SelectSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[selectfile,selectpath]=uigetfile('*.mat','Select Settings file');
if selectfile~=0
    set(handles.ET_load_path,'string',[selectpath '\' selectfile]);
end

% --- Executes on button press in SelectRestart.
function SelectRestart_Callback(hObject, eventdata, handles)
% hObject    handle to SelectRestart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[selectfile,selectpath]=uigetfile('Select Restart files','MultiSelect','on');
if selectfile==0
    return
end
if ~iscell(selectfile)
    restartstring=[selectpath selectfile];
else
    restartstring=[];
    if strcmp(selectpath,[pwd '\'])
        selectpath='';
    end
    for i_f=1:size(selectfile,2)-1
       restartstring=[restartstring selectpath selectfile{1,i_f} ';'];
    end
    restartstring=[restartstring selectpath selectfile{1,end}]; 
end
set(handles.ET_path_for_restart,'string',restartstring);

% --- Executes on button press in SelectHier.
function SelectHier_Callback(hObject, eventdata, handles)
% hObject    handle to SelectHier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[selectfile,selectpath]=uigetfile('*.mat','Select Hier file');
set(handles.ET_hier_model_mat,'string',[selectpath '\' selectfile]);

function UpdateFiles(n_ch_out,DoIndividualFiles)
%EFF21: Update seizure-info (do this once per subfile). THESE ARE SPARSE UPDATES
%EFF21, Also, faster to use save with append if you are not growing dynamically
global mf mf2 save_mat_path save_mat_path2
global Seizure_Count Seizure_Count_Current n_read
global Sz_On_Index Sz_Off_Index Sz_Index
global Tr_On_Index Tr_Off_Index Tr_Index
global spike_count_history spike_count_history_raw spike_count_history_combined 
global ro_data ro_dd ro_ard  %indices of final rows in data, detect_data and art_rem_data respectively
global DoTroubleShoot
global AlreadyUpdated

AlreadyUpdated=1; %Change flag so that we are not updating twice for same file (relevant for DoExisting (aka offline ) mode)

if DoTroubleShoot
    disp(['Run: UpdateFiles (' num2str(DoIndividualFiles) ')...'])    
end

%save to summary file and individual files
if ~isempty(Sz_On_Index)
    stim_flag=sparse(Sz_On_Index(:,1),Sz_On_Index(:,2),true(size(Sz_On_Index,1),1));
    if DoIndividualFiles
        Sz_On=sparse(Sz_On_Index(Sz_Index(1):end,1),Sz_On_Index(Sz_Index(1):end,2),true(size(Sz_On_Index,1)-Sz_Index(1)+1,1));
    end
else
    stim_flag=sparse([]);
    Sz_On=sparse([]);
end
if ~isempty(Sz_Off_Index)
    Seizure_Off=sparse(Sz_Off_Index(:,1),Sz_Off_Index(:,2),true(size(Sz_Off_Index,1),1));
    if DoIndividualFiles
        Sz_Off=sparse(Sz_Off_Index(Sz_Index(2):end,1),Sz_Off_Index(Sz_Index(2):end,2),true(size(Sz_Off_Index,1)-Sz_Index(2)+1,1));
    end
else
    Seizure_Off=sparse([]);
    Sz_Off=sparse([]);
end
if ~isempty(Tr_On_Index)
    spikes_trigger=sparse(Tr_On_Index(:,1),Tr_On_Index(:,2),true(size(Tr_On_Index,1),1));
    if DoIndividualFiles
        Tr_On=sparse(Tr_On_Index(Tr_Index(1):end,1),Tr_On_Index(Tr_Index(1):end,2),true(size(Tr_On_Index,1)-Tr_Index(1)+1,1));
    end
else
    spikes_trigger=sparse([]);
    Tr_On=sparse([]);
end
if ~isempty(Tr_Off_Index)
    spikes_trigger_stop=sparse(Tr_Off_Index(:,1),Tr_Off_Index(:,2),true(size(Tr_Off_Index,1),1));
    if DoIndividualFiles
        Tr_Off=sparse(Tr_Off_Index(Tr_Index(2):end,1),Tr_Off_Index(Tr_Index(2):end,2),true(size(Tr_Off_Index,1)-Tr_Index(2)+1,1));
    end
else
    spikes_trigger_stop=sparse([]);
    Tr_Off=sparse([]);
end

%Update bookkeeping indices for sparse arrays
Sz_Index(1)=size(Sz_On_Index,1)+1;      %Index where seizure info for current file will start
Sz_Index(2)=size(Sz_Off_Index,1)+1;     %Index where seizure info for current file will start
Tr_Index(1)=size(Tr_On_Index,1)+1;
Tr_Index(2)=size(Tr_Off_Index,1)+1;    

%save updated variables in one call to save with append (faster than using matfile. repeatedly
save(save_mat_path,'Seizure_Count','n_ch_out','spikes_trigger','spikes_trigger_stop','stim_flag','Seizure_Off','-append');
if DoIndividualFiles
    save(save_mat_path2,'Seizure_Count_Current','Tr_On','Tr_Off','Sz_On','Sz_Off','ro_data','ro_dd','ro_ard','-append');
end

%save final spike count since last save, this will be appended, since it is
%part of bigger structure
if mod(n_read,10)~=0 
    mf.spike_count(n_read-mod(n_read,10)+1:n_read,:) = [spike_count_history(11-mod(n_read,10):10,:) ...          %spike_count_history
                                                                spike_count_history_raw(11-mod(n_read,10):10,:) ...      %spike_count_raw
                                                                spike_count_history_combined(11-mod(n_read,10):10,:)];   %spike_count_combined
end
    

% --- Executes on button press in SelectReanalysisDir.
function SelectReanalysisDir_Callback(hObject, eventdata, handles)
% hObject    handle to SelectReanalysisDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ipd=get(handles.InputFileName,'string');
if exist(ipd,'dir')
    selectpath = uigetdir(ipd,'Select location for reanalysis files');
else
    selectpath = uigetdir(pwd,'Select location for reanalysis files');
end
if selectpath~=0
    set(handles.InputFileName,'string',selectpath);
end
