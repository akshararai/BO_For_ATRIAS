function nm = fn_set_vnmc_params(sm, sample_time)


% 
% ATRIAS_VMSMlInit.m
%   set parameters of the virtual musculoskeletal model
%
% set:
%   segment dimensions
%   muscle-skeleton attachments
%   muscle dynamics
%   neural transmission delays
%
% by Seungmoon Song
% started at Oct 28, 2014
%


% ================= %
% VIRTUAL-LEG MODEL %
% ================= %

% segment lengths
nm.l_seg = sm.l_seg;
nm.min_splay_angle = 45*pi/180; % DEFINES STRAIGHT LEG LENGTH delta_phi_min = 25.5*pi/180;
nm.l_v = sm.l_seg*cos(nm.min_splay_angle/2);  %CS1 (ankle joint) to C2 (knee joint)

% =========================== %
% MUSCLE-SKELETON ATTACHMENTS %
% =========================== %
       
% Hip FLexor group attachement
nm.rHFL       =       0.08; % [m]   constant lever contribution 
nm.phirefHFL  = 160*pi/180; % [rad] reference angle at which MTU length equals 
nm.rhoHFL     =        0.5; %       sum of lopt and lslack          

% GLUtei group attachement
nm.rGLU       =       0.08; % [m]   constant lever contribution 
nm.phirefGLU  = 120*pi/180; % [rad] reference angle at which MTU length equals 
nm.rhoGLU     =        0.5; %       sum of lopt and lslack 
                         
% HAMstring group attachement (hip)
nm.rHAMh       = 0.08;         % [m]   constant lever contribution 
nm.phirefHAMh  = 150*pi/180;   % [rad] reference angle at which MTU length equals 
nm.rhoHAMh     = 0.5;          %       sum of lopt and lslack 

% HAMstring group attachement (knee)
nm.rHAMk       = 0.05;         % [m]   constant lever contribution 
nm.phirefHAMk  = 180*pi/180;   % [rad] reference angle at which MTU length equals 
nm.rhoHAMk     = 0.5;          %       sum of lopt and lslack 

% RF group attachement (hip)
nm.rRFh      =       0.08; % [m]   constant lever contribution 
nm.phirefRFh = 170*pi/180; % [rad] reference angle at which MTU length equals 
nm.rhoRFh    =        0.3; %       sum of lopt and lslack 

% RF group attachement (knee)
nm.rRFkmax     = 0.06;         % [m]   maximum lever contribution
nm.rRFkmin     = 0.04;         % [m]   minimum lever contribution
nm.phimaxRFk   = 165*pi/180;   % [rad] angle of maximum lever contribution
nm.phiminRFk   =  45*pi/180;   % [rad] angle of minimum lever contribution
nm.phirefRFk   = 125*pi/180;   % [rad] reference angle at which MTU length equals 
nm.rhoRFk      = 0.5;          %       sum of lopt and lslack 
nm.phiScaleRFk = acos(nm.rRFkmin/nm.rRFkmax)/(nm.phiminRFk-nm.phimaxRFk);

% VAStus group attachement
nm.rVASmax     = 0.06;         % [m]   maximum lever contribution
nm.rVASmin     = 0.04;         % [m]   minimum lever contribution
nm.phimaxVAS   = 165*pi/180;   % [rad] angle of maximum lever contribution
nm.phiminVAS   =  45*pi/180;   % [rad] angle of minimum lever contribution
nm.phirefVAS   = 120*pi/180;   % [rad] reference angle at which MTU length equals 
nm.rhoVAS      = 0.6;          %       sum of lopt and lslack
nm.phiScaleVAS = acos(nm.rVASmin/nm.rVASmax)/(nm.phiminVAS-nm.phimaxVAS);

% BFSH group attachement
nm.rBFSH    	= 0.04;         % [m]   constant lever contribution 
nm.phirefBFSH 	= 160*pi/180;   % [rad] reference angle at which MTU length equals 
nm.rhoBFSH    	= 0.7;          %       sum of lopt and lslack


% =============== %
% MUSCLE DYNAMICS %
% =============== %

% -------------------------------
% shared muscle tendon parameters
% -------------------------------

% excitation-contraction coupling
nm.preA =  0.01; %[] preactivation
nm.tau  =  0.01; %[s] delay time constant

% contractile element (CE) force-length relationship
nm.w    =   0.56; %[lopt] width
nm.c    =   0.05; %[]; remaining force at +/- width

% CE force-velocity relationship
nm.N    =   1.5; %[Fmax] eccentric force enhancement
nm.K    =     5; %[] shape factor

% Series elastic element (SE) force-length relationship
nm.eref =  0.04; %[lslack] tendon reference strain

% -------------------------------
% shared muscle tendon parameters (DISCRETE)
% -------------------------------
% K_ECC_discrete = sample_time/((2*tau/sample_time)+sample_time);
% a_ECC_discrete = -1;
% b_ECC_discrete = ((2*tau/sample_time)-sample_time)/((2*tau/sample_time)+sample_time);


% --------------------------
% Soft Joint Limits
% --------------------------
nm.v12_max = 179*pi/180; %max virtual knee angle (to avoid singularities)

% angles at which soft limits engages
nm.phi12_min = 0*pi/180; %[rad]
nm.phi12_max = 175*pi/180; %[rad]
    
nm.phi23_max = 230*pi/180; %[rad]
nm.phi23_min = 100*pi/180; %[rad]

% soft block reference joint stiffness
nm.c_jointstop     = 0.3 / (pi/180);  %[Nm/rad]

% soft block maximum joint stop relaxation speed
nm.w_max_jointstop = 1 * pi/180; %[rad/s]

% --------------------------
% muscle-specific parameters
% --------------------------

% hip flexor muscles
nm.FmaxHFL   = 2000; % maximum isometric force [N]
nm.loptHFL   = 0.11; % optimum fiber length CE [m]
nm.vmaxHFL   =   12; % maximum contraction velocity [lopt/s]
nm.lslackHFL = 0.10; % tendon slack length [m]

% glutei muscles
nm.FmaxGLU   = 1500; % maximum isometric force [N]
nm.loptGLU   = 0.11; % optimum fiber length CE [m]
nm.vmaxGLU   =   12; % maximum contraction velocity [lopt/s]
nm.lslackGLU = 0.13; % tendon slack length [m]

% hamstring muscles
nm.FmaxHAM   = 3000; % maximum isometric force [N]
nm.loptHAM   = 0.10; % optimum fiber length CE [m]
nm.vmaxHAM   =   12; % maximum contraction velocity [lopt/s]
nm.lslackHAM = 0.31; % tendon slack length [m]

% rectus femoris muscles
nm.FmaxRF   = 1200; % %850 maximum isometric force [N]
nm.loptRF   = 0.08; % optimum fiber length CE [m]
nm.vmaxRF   =   12; % maximum contraction velocity [lopt/s]
nm.lslackRF = 0.35; % tendon slack length [m]

% vasti muscles
nm.FmaxVAS     = 6000; % maximum isometric force [N]
nm.loptVAS     = 0.08; % optimum fiber length CE [m]
nm.vmaxVAS     =   12; % maximum contraction velocity [lopt/s]
nm.lslackVAS   = 0.23; % tendon slack length [m]

% BFSH
nm.FmaxBFSH	=  350; % maximum isometric force [N]
nm.loptBFSH    = 0.12; % optimum fiber length CE [m]
nm.vmaxBFSH    =   12; %6 % maximum contraction velocity [lopt/s]
nm.lslackBFSH  = 0.10; % tendon slack length [m]

% gastrocnemius muscle
nm.FmaxGAS    = 1500; % maximum isometric force [N]
nm.loptGAS    = 0.05; % optimum fiber length CE [m]
nm.vmaxGAS    =   12; % maximum contraction velocity [lopt/s]
nm.lslackGAS  = 0.40; % tendon slack length [m]

% soleus muscle
nm.FmaxSOL    = 4000; % maximum isometric force [N]
nm.loptSOL    = 0.04; % optimum fiber length CE [m]
nm.vmaxSOL    =    6; % maximum contraction velocity [lopt/s]
nm.lslackSOL  = 0.26; % tendon slack length [m]

% tibialis anterior
nm.FmaxTA     =  800; % maximum isometric force [N]
nm.loptTA     = 0.06; % optimum fiber length CE [m]
nm.vmaxTA     =   12; % maximum contraction velocity [lopt/s]
nm.lslackTA   = 0.24; % tendon slack length [m]


% --------------------------------------------------
% 5 Muscle parameters for Calculating Metabolic Cost
% --------------------------------------------------
% added by S. Song, May 2011

nm.FTproportionHFL = 0.45;
nm.massHFL = 0.95; % kg

nm.FTproportionGLU = 0.50;
nm.massGLU = 0.70; % kg

nm.FTproportionHAM = 0.40;
nm.massHAM = 1.25; % kg

nm.FTproportionRF = 0.55;
nm.massRF = 0.40; % kg

nm.FTproportionVAS = 0.50;
nm.massVAS = 2.05; % kg

nm.FTproportionBFSH = 0.35;
nm.massBFSH = 0.15; % kg

% FTproportionGAS = 0.50;
% massGAS = 0.30; % kg
% 
% FTproportionTA = 0.30;
% massTA = 0.20; % kg
% 
% FTproportionSOL = 0.25;
% massSOL = 0.70; % kg



% ========================== %
% NEURAL TRAJSMISSION DELAYS %
% ========================== %

% nm.LongLoopDelay	= 0.030;    % [s] additional to spinal reflexes
% nm.LongDelay       = 0.020/2;  % [s] ankle joint muscles
nm.MidDelay        = 0.010/2;  % [s] knee joint muscles
% nm.ShortDelay      = 0.005/2;  % [s] hip joint muscles
% nm.MinDelay        = 0.001/2;  % [s] between neurons in the spinal cord

% for descrete time step controller

% f = round(LongLoopDelay/sample_time);
% LongDelay_discrete = round(LongDelay/sample_time);
nm.MidDelay_discrete = round(nm.MidDelay/sample_time);
% ShortDelay_discrete = round(ShortDelay/sample_time);
% MinDelay_discrete = round(MinDelay/sample_time);



