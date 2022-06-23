% Get 4D references w/ search-space constrained to real-data
% characteristics
%% Init
data_dir = 'data\measured_data\awa_pm_45\cT_4D.mat';

if isfile(data_dir)
    load(data_dir);
else
    disp('Data_dir is not file.\n')
end

% [WIP]
% Get local copies
cT             = data.cT;
sheeting_angle = data.sheeting_angle;
AWA            = data.AWA;

L_sa    = size(data.sheeting_angle, 2);
L_awa   = size(data.AWA, 2);

% Concatenate axes
AWA_cat = [-flip(AWA), AWA];

sheeting_angle_cat                    = zeros(4, L_sa, L_awa*2);
sheeting_angle_cat(:, :, 1:L_awa)     = -flip(sheeting_angle, 3); % flip in AWA
sheeting_angle_cat(:, :, L_awa+1:end) = sheeting_angle;

% Concatenate cT
cT_cat                         = zeros(L_awa*2, L_sa, L_sa, L_sa, L_sa);
cT_cat(1:L_awa, :, :, :, :)    = flip(cT); % flip in = AWA
cT_cat(L_awa+1:end, :, :, :, :) = cT;
