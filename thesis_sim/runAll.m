%% RUN ALL THESIS EXPERIMENTS (H1, H2, sensitivity S)
% Results in results/, figures in figures/. Each experiment script has a
% "rerun" flag at the top: set it to false to only redraw the figures.

exp_H1;     % radius-scheduled mission speed vs original Durmaz2024
exp_H2;     % safety filters under unmodeled disturbances
exp_S;      % sensitivity of the radius-scheduled mission speed
