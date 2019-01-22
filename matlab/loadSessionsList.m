function [columnLabels S] = loadSessionsList
%function [sessions doses group] = loadSessionsList
%
% Function that loads the details about the experiment conducted by Diego
% H. de Miranda during his MSc. project.
% UFABC - 2014
% Each drug was injected one or more times. Sessions with the same name
% were replicaton of the experiment for precision increase
% Marcelo Bussotti Reyes - UFABC 2014


% halo_low  = haloperidol 0.2 mg/kg
% halo_high = haloperidol 1.0 mg/kg (two sessions)
% olan_high = 
%


% To-do list:
% 1.Check if colHalo1mg = 6 bellow is correct. It seems that the low dose 
%   is 0.2mg/kg
% 2.[solved, needs implementation] check if it is straighforward to load
%   from excel sheet.
%   use xlsread, but there must be a simple worksheet, no comments no
%   nothing. There may be a header with lables. To load this, use [header
%   data] = xlsread (check the other options). It did not work for xlsx
%   sheets.
%3. Define what to do when there are two sessions for the same dose.
%   Options: merge, or define the session to be used.
%__________________________________________________________________________

% Code for experiment identification, which is also the prefix for the data
% file names
%prefix = 'AB1';

columnLabels = {'animal_id'	'baseline'	'baseline_2'	'baseline_ic' 'saline_ip'...        % 1 - 5
    ''              'saline_ic'     'etanol_ic'     ''              'saline+tween_ic'...    % 6 - 10
    'halo_4_ic'     'halo_4_ic'     'halo_1_ic'     'halo_5_ic'     'halo_20_ic'...         % 11 - 15
    'halo_1_ip'     ''          	'halo_0.2_ip'	'halo_0.07_ip'	'olan_3_ip'...          % 16 - 20
    ''          	'olan_2_ip'     'olan_1_ip'     'olan_0.3_ip'	'olan_1_ic'...          % 21 - 25
    'olan_5_ic'     'olan_10_ic'	'apo_1_ip'      'apo_3_ip'      'apo_9_ip'...           % 26 - 30
    'apo_0.3_ip'	'apo_1_ic'      'apo_5_ic'      'apo_10_ic'     'amph_2_ic'};           % 31 - 35

% session matrix S: contains all the subjects and sessions for all subjects
% S = [...
% 7	17	21	8	18	31	6	0	0	26	27	9	7	0	0	19	23	29	0	37	0	35	33	0	0	0	0	0	0	0	0	0	0	0	15
% 8	17	21	8	18	31	6	0	0	0	0	9	7	0	0	19	23	29	0	37	0	35	33	0	0	0	0	0	0	0	0	0	0	0	15
% 9	17	21	8	18	31	6	0	0	26	27	9	7	0	0	19	23	29	0	37	0	35	33	0	0	0	0	0	0	0	0	0	0	0	15
% 10	17	21	8	18	31	6	0	0	26	27	9	7	0	0	19	23	29	0	37	0	35	33	0	0	0	0	0	0	0	0	0	0	0	15
% 11	17	21	8	18	31	6	0	0	0	0	9	7	0	0	19	23	29	0	37	0	35	33	0	0	0	0	0	0	0	0	0	0	0	15
% 12	17	21	8	18	31	6	0	0	26	27	9	7	0	0	19	23	29	0	37	0	35	33	0	0	0	0	0	0	0	0	0	0	0	15
% 13	17	21	8	18	31	6	0	0	26	27	9	7	0	0	19	23	29	0	37	0	35	33	0	0	0	0	0	0	0	0	0	0	0	15
% 14	17	21	8	18	31	6	0	0	26	27	9	7	0	0	19	23	29	0	37	0	35	33	0	0	0	0	0	0	0	0	0	0	0	15
% 15	17	21	8	18	31	6	0	0	26	27	9	7	0	0	19	23	29	0	37	0	35	33	0	0	0	0	0	0	0	0	0	0	0	15
% 34	17	27	0	18	28	32	0	0	0	0	0	0	0	0	20	24	22	26	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
% 35	17	27	0	18	28	0	0	0	0	0	0	0	0	0	22	24	20	26	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
% 36	17	0	44	18	0	32	42	47	0	0	0	35	37	39	20	24	22	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	45	0
% 37	17	27	0	18	28	0	0	0	0	0	0	0	0	0	0	0	0	0	22	24	0	20	26	0	0	0	0	0	0	0	0	0	0	0
% 38	17	27	0	18	28	0	0	0	0	0	0	0	0	0	0	0	0	0	20	24	0	22	26	0	0	0	0	0	0	0	0	0	0	0
% 39	17	27	0	18	28	0	0	0	0	0	0	0	0	0	0	0	0	0	22	24	0	20	26	0	0	0	0	0	0	0	0	0	0	0
% 40	17	27	0	18	28	0	0	0	0	0	0	0	0	0	0	0	0	0	20	24	0	22	26	0	0	0	0	0	0	0	0	0	0	0
% 41	17	27	0	18	28	0	0	0	0	0	0	0	0	0	0	0	0	0	22	24	0	20	26	0	0	0	0	0	0	0	0	0	0	0
% 42	17	27	44	18	28	32	42	50	0	0	0	0	0	45	0	0	0	0	20	24	0	22	26	35	37	39	0	0	0	0	0	0	47	0
% 43	17	27	0	18	28	32	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	22	20	24	26	0	0	0	0
% 44	17	27	0	18	28	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	20	22	24	26	0	0	0	0
% 45	17	27	44	18	28	32	42	47	0	0	0	0	0	45	0	0	0	0	0	0	0	0	0	0	0	0	22	20	24	26	35	37	39	0
% 46	17	27	0	18	28	0	0	0	0	0	0	0	0	0	22	24	20	26	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
% 47	17	27	0	18	28	0	0	0	0	0	0	0	0	0	20	24	22	26	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
% 48	17	27	44	18	28	32	42	47	0	0	0	35	37	39	22	24	20	26	0	0	0	0	0	0	0	0	0	0	0	0	0	0	45	0
% 49	17	27	0	18	28	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	20	22	24	26	0	0	0	0
% 50	17	27	44	18	28	32	42	47	0	0	0	0	0	45	0	0	0	0	0	0	0	0	0	0	0	0	22	20	24	26	35	37	39	0
% 51	17	27	0	18	28	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	20	22	24	26	0	0	0	0
% 52	17	27	44	18	28	32	42	47	0	0	0	35	37	39	0	0	0	0	20	24	0	22	26	0	0	0	0	0	0	0	0	0	45	0
% 53	17	27	44	18	28	32	42	47	0	0	0	35	37	39	0	0	0	0	20	24	0	22	26	0	0	0	0	0	0	0	0	0	45	0
% 54	17	27	44	18	28	32	42	47	0	0	0	35	37	39	0	0	0	0	22	24	0	20	26	0	0	0	0	0	0	0	0	0	45	0
% 55	17	27	44	18	28	32	42	50	0	0	0	0	0	47	20	24	22	26	0	0	0	0	0	35	37	39	0	0	0	0	0	0	45	0
% 56	17	27	44	18	28	32	42	50	0	0	0	0	0	47	20	24	22	26	0	0	0	0	0	35	37	39	0	0	0	0	0	0	45	0
% 57	17	27	44	18	28	32	42	50	0	0	0	0	0	47	22	24	20	26	0	0	0	0	0	35	37	39	0	0	0	0	0	0	45	0
% 58	17	27	0	18	28	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	22	20	24	26	0	0	0	0
% 59	17	27	44	18	28	32	42	47	0	0	0	0	0	45	0	0	0	0	0	0	0	0	0	0	0	0	22	20	24	26	35	37	39	0
% 60	17	27	44	18	28	32	42	47	0	0	0	0	0	45	0	0	0	0	0	0	0	0	0	0	0	0	20	22	24	26	35	37	39	0
% 61	17	27	44	18	28	32	42	50	0	0	0	0	0	45	0	0	0	0	22	24	0	20	26	35	37	39	0	0	0	0	0	0	47	0
% 62	17	27	44	18	28	32	42	47	0	0	0	0	0	45	0	0	0	0	0	0	0	0	0	0	0	0	20	22	24	26	35	37	39	0
% 63	17	27	44	18	28	32	42	50	0	0	0	0	0	45	22	24	20	26	0	0	0	0	0	35	37	39	0	0	0	0	0	0	47	0
% ];

S = [...
34	17	26	0	18	28	32	0	0	0	0	0	0	0	0	20	24	22	26	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
35	17	26	0	18	28	0	0	0	0	0	0	0	0	0	22	24	20	26	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
36	17	0	38	18	0	32	42	47	0	0	0	35	37	39	20	24	22	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	45	0
37	17	26	0	18	28	0	0	0	0	0	0	0	0	0	0	0	0	0	22	24	0	20	26	0	0	0	0	0	0	0	0	0	0	0
38	17	26	0	18	28	0	0	0	0	0	0	0	0	0	0	0	0	0	20	24	0	22	26	0	0	0	0	0	0	0	0	0	0	0
39	17	26	0	18	28	0	0	0	0	0	0	0	0	0	0	0	0	0	22	24	0	20	26	0	0	0	0	0	0	0	0	0	0	0
40	17	26	0	18	28	0	0	0	0	0	0	0	0	0	0	0	0	0	20	24	0	22	26	0	0	0	0	0	0	0	0	0	0	0
41	17	26	0	18	28	0	0	0	0	0	0	0	0	0	0	0	0	0	22	24	0	20	26	0	0	0	0	0	0	0	0	0	0	0
42	17	26	38	18	28	32	42	50	0	0	0	0	0	45	0	0	0	0	20	24	0	22	26	35	37	39	0	0	0	0	0	0	47	0
43	17	26	0	18	28	32	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	22	20	24	26	0	0	0	0
44	17	26	0	18	28	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	20	22	24	26	0	0	0	0
45	17	26	38	18	28	32	42	47	0	0	0	0	0	45	0	0	0	0	0	0	0	0	0	0	0	0	22	20	24	26	35	37	39	0
46	17	26	0	18	28	0	0	0	0	0	0	0	0	0	22	24	20	26	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
47	17	26	0	18	28	0	0	0	0	0	0	0	0	0	20	24	22	26	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
48	17	26	38	18	28	32	42	47	0	0	0	35	37	39	22	24	20	26	0	0	0	0	0	0	0	0	0	0	0	0	0	0	45	0
49	17	26	0	18	28	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	20	22	24	26	0	0	0	0
50	17	26	38	18	28	32	42	47	0	0	0	0	0	45	0	0	0	0	0	0	0	0	0	0	0	0	22	20	24	26	35	37	39	0
51	17	26	0	18	28	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	20	22	24	26	0	0	0	0
52	17	26	38	18	28	32	42	47	0	0	0	35	37	39	0	0	0	0	20	24	0	22	26	0	0	0	0	0	0	0	0	0	45	0
53	17	26	38	18	28	32	42	47	0	0	0	35	37	39	0	0	0	0	20	24	0	22	26	0	0	0	0	0	0	0	0	0	45	0
54	17	26	38	18	28	32	42	47	0	0	0	35	37	39	0	0	0	0	22	24	0	20	26	0	0	0	0	0	0	0	0	0	45	0
55	17	26	38	18	28	32	42	50	0	0	0	0	0	47	20	24	22	26	0	0	0	0	0	35	37	39	0	0	0	0	0	0	45	0
56	17	26	38	18	28	32	42	50	0	0	0	0	0	47	20	24	22	26	0	0	0	0	0	35	37	39	0	0	0	0	0	0	45	0
57	17	26	38	18	28	32	42	50	0	0	0	0	0	47	22	24	20	26	0	0	0	0	0	35	37	39	0	0	0	0	0	0	45	0
58	17	26	0	18	28	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	22	20	24	26	0	0	0	0
59	17	26	38	18	28	32	42	47	0	0	0	0	0	45	0	0	0	0	0	0	0	0	0	0	0	0	22	20	24	26	35	37	39	0
60	17	26	38	18	28	32	42	47	0	0	0	0	0	45	0	0	0	0	0	0	0	0	0	0	0	0	20	22	24	26	35	37	39	0
61	17	26	38	18	28	32	42	50	0	0	0	0	0	45	0	0	0	0	22	24	0	20	26	35	37	39	0	0	0	0	0	0	47	0
62	17	26	38	18	28	32	42	47	0	0	0	0	0	45	0	0	0	0	0	0	0	0	0	0	0	0	20	22	24	26	35	37	39	0
63	17	26	38	18	28	32	42	50	0	0	0	0	0	45	22	24	20	26	0	0	0	0	0	35	37	39	0	0	0	0	0	0	47	0
];

% % --- Groups organized by the drugs that were injected ---
% halo = 1; 
% olan = 2;
% apo  = 3;
% sessions = {'halo' 'olan' 'apo'};
% 
% % --- doses --- 
% base = 1;
% sal  = 2;
% d1mg = 3;       % 1mg/kg (low)
% d3mg = 4;       % 3mg/kg (high)
% d9mg = 5;       % 9mg/kg (super)
% 
% doses    = {'base' 'sal'  '1mg' '3mg' '9mg'};
% 
% % --- Animals that belong to each group ---
% % this routine looks for the animals (column 1) that had sessions 
% % (different from zero) and stores them in the subjs sub-structure
% 
% % --- Haloperidol ---
% group(halo,d1mg).subjs = S(S(:,colHalo1mg)>0, colSubjs);
% group(halo,d3mg).subjs = S(S(:,colHalo3mg)>0, colSubjs);
% if ~(group(halo,d1mg).subjs==group(halo,d3mg).subjs)
% 	disp('Attention: halo 1 e 3mg: different number of subjcts');
% end
% group(halo,base).subjs = S(S(:,colHalo1mg)>0, colBaseline);
% group(halo,sal ).subjs = S(S(:,colHalo1mg)>0, colSaline);
% 
% % --- Olanzapina ---
% group(olan,d1mg).subjs = S(S(:,colOlan1mg)>0, colSubjs);
% group(olan,d3mg).subjs = S(S(:,colOlan3mg)>0, colSubjs);
% if ~(group(olan,d1mg).subjs==group(olan,d3mg).subjs)
% 	disp('Attention: olan 1 e 3mg: Número de sujeitos diferente');
% end
% group(olan,base).subjs = S(S(:,colOlan1mg)>0, colBaseline);
% group(olan,sal ).subjs = S(S(:,colOlan1mg)>0, colSaline);
% 
% % --- apomorfina ---
% group(apo, d1mg).subjs = S(S(:,colApo1mg)>0, colSubjs);
% group(apo, d3mg).subjs = S(S(:,colApo3mg)>0, colSubjs);
% if group(apo,d1mg).subjs ~= group(apo,d3mg).subjs
% 	disp('Attention: apomorphine 1 and 3mg: mumber of subjects is diferent');
% end
% group(apo, base).subjs = S(S(:,colApo1mg)>0, colBaseline);
% group(apo, sal ).subjs = S(S(:,colApo1mg)>0, colSaline);
% 
% % --- Now the drug sessions (might be more than one session for the same
% % dose)
% group(halo,d1mg).sessions = S(S(:,colHalo1mg)>0, colHalo1mg);
% group(halo,d3mg).sessions = S(S(:,colHalo3mg)>0, colHalo3mg);
% 
% %group(halo).dose(d3mg).S = S(S(:,colHalo3mg)>0, colHalo3mg);
% 
% 
% 
% 
% 
% 
% 
