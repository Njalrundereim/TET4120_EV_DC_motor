function [sys,x0,str,ts] = sfun_onephase_2LDC_A(t,x,u,flag,theta_phase,t_dead,t_min,type_signal,asym)
%   Variable step S-function.
%   Simulation of one leg of the pwm-module, including non-linear
%   effects as blanking time and minimum turn-off time
switch flag,
   
  case 0,    % Initialization %
    [sys,x0,str,ts]=mdlInitializeSizes;
  
  case 3,    % Outputs %
    sys=mdlOutputs(t,x,u);

  case 4,    % GetTimeOfNextVarHit %
    sys=mdlGetTimeOfNextVarHit(t,x,u,theta_phase,t_dead,t_min,type_signal,asym) ;
  case 9,    % Terminate %
    sys=mdlTerminate(t,x,u); 
  
  case {1,2}    % Unhandled flags %
    sys = []; 
  
  otherwise    % Unexpected flags %
    error(['Unhandled flag = ',num2str(flag)]);

end

% end sfun_onephase

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts]=mdlInitializeSizes

sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 5;    
sizes.NumInputs      = 4;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;     % at least one sample time is needed

sys = simsizes(sizes);
clc
%
% initialize the initial conditions
%
x0  = [];
str = [];
k = -1;    			  % starting with -1 gives k=0 in the first step
dt_plot = 1.0e-8;   % sets the time between two following calculation steps
						  % needed to get high quality plots of the digital signals
%                     
% Initialize the UserData field of the block
%
set_param(gcbh,'UserData',[k,dt_plot,0     ,0     ,0     ,0       ,0        ,0      ,0       ,0   ,0     ,0      ,0        ,0]);
% vector is defined by: % [k,dt_plot,t_sync,t_next,t_xtra,t_a_u_on,t_a_u_off,t_a_lon,t_a_loff,u_st_a,status,status_,dt_a_u_on,ignore_xtra]
%
%
% initialize the array of sample times
%
ts  = [-2 0];      % variable sample time

% end mdlInitializeSizes

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u);

% the if-statement defines the s-function's digital output  
% dependent of the state of the counter vaiable k

variable_vect = get_param(gcbh,'UserData');
k = variable_vect(1);
dt_a_u_on = variable_vect(13);
t_sync = variable_vect(3);
status = variable_vect(11); %Saturation status upper band. 1 = overmodulation, 0 = undermodulation, 2 = no saturation
status_ = variable_vect(12); %Previous saturation status upper band
ignore_xtra = variable_vect(14);
sys(5) = variable_vect(10);
if k == 0 %Start of trigperiod
  if ((status == 1 && status_ == 1) || (status_==1 && ignore_xtra ==1)) %Positive --> Positive saturation
      sys(1)=1;
      sys(2)=0;     
  elseif status_==1 %Positive --> Negative or no saturation
      sys(1)=0;
      sys(2)=0;
  else  
      sys(1)=0;
      sys(2)=1;
  end
  sys(3) = t_sync;% ------------- '' ------------- , t_sync is unchanched
  sys(4) = 0;     % ------------- '' ------------- , sync is (still) equal to 0
elseif k == 1 %Extra switching instant - dt_plot (t_xtra-dt_plot)
  if ((status == 1 && status_ == 1) || (status_==1 && ignore_xtra ==1))
      sys(1)=1;
      sys(2)=0 ;  
  elseif (status_==1)  %Positive --> Negative or no saturation
      sys(1)=0;
      sys(2)=0;
  else
      sys(1)=0;
      sys(2)=1;
  end
  sys(3) = t_sync;% ------------- '' ------------- , t_sync is unchanched
  sys(4) = 0;     % ------------- '' ------------- , sync is (still) equal to 0
elseif k == 2 %Extra switching instant (t_xtra)
  if ((status == 1 && status_ == 1) || (status_==1 && ignore_xtra ==1)) %Positive to positive saturation
      sys(1)=1;
      sys(2)=0 ;  
  else
      sys(1)=0;
      sys(2)=1;
  end
  sys(3) = t_sync;% ------------- '' ------------- , t_sync is unchanched
  sys(4) = 0;     % ------------- '' ------------- , sync is (still) equal to 0
elseif k == 3 %t_a_l_off - dt_plot
  if ((status == 1 && status_ == 1) || (status_==1 && ignore_xtra ==1)) %Positive saturation --> Positive Saturation
      sys(1)=1;
      sys(2)=0;
  else
      sys(1)=0;
      sys(2)=1;
  end
  sys(3) = t_sync;% ------------- '' ------------- , t_sync is unchanched
  sys(4) = 0;     % ------------- '' ------------- , sync is (still) equal to 0
elseif k == 4 %t_a_l_off
  if ((status == 1 && status_ == 1) || (status_==1 && ignore_xtra ==1)) %Positive saturation --> Positive Saturation
      sys(1)=1;
      sys(2)=0;
  elseif status == 0 %  --> Negative Saturation
      sys(1)=0;
      sys(2)=1;
  else
      sys(1)=0; 
      sys(2)=0;
  end
  sys(3) = t_sync;% ------------- '' ------------- , t_sync is unchanched
  sys(4) = 0;     % ------------- '' ------------- , sync is (still) equal to 0
elseif k == 5 %t_a_u_on - dt_plot
  if ((status == 1 && status_ == 1) || (status_==1 && ignore_xtra ==1)) %Positive saturation --> Positive Saturation
      sys(1)=1;
      sys(2)=0;
  elseif status == 0 %  --> Negative Saturation
      sys(1)=0;
      sys(2)=1;
  else
      sys(1)=0; 
      sys(2)=0;
  end
  sys(3) = t_sync;% ------------- '' ------------- , t_sync is unchanched
  sys(4) = 0;     % ------------- '' ------------- , sync is (still) equal to 0
elseif k == 6 %t_a_u_on
  if status == 0  % -->Negative Saturation
      sys(1)=0;
      sys(2)=1;
  else
      sys(1)=1;
      sys(2)=0;
  end
  sys(3) = t_sync;% ------------- '' ------------- , t_sync is unchanched
  sys(4) = 0;     % ------------- '' ------------- , sync is (still) equal to 0
elseif k == 7 %t_sync - dt_plot
  if status == 0  % -->Negative Saturation
      sys(1)=0;
      sys(2)=1;
  else
      sys(1)=1;
      sys(2)=0;
  end
  sys(3) = t_sync;% ------------- '' ------------- , t_sync is unchanched
  sys(4) = 0;     % ------------- '' ------------- , sync is (still) equal to 0
  
elseif k == 8 %t_sync
  if status == 0  % -->Negative Saturation
      sys(1)=0;
      sys(2)=1;
  else
      sys(1)=1;
      sys(2)=0;
  end
  sys(3) = t_sync;% --------- '' --------- , t_sync is unchanched
  sys(4) = 1;     % --------- '' --------- , sync is set equal to 1
  %sys(5) = abs(-(4/u(3))*(t-round(t/u(3))*u(3)))-1;
  
elseif k == 9 %t_a_u_off-dt_plot
  if status == 0  % -->Negative Saturation
      sys(1)=0;
      sys(2)=1;
  else
      sys(1)=1;
      sys(2)=0;
  end
  sys(3) = t_sync;% -------------- '' --------------- , t_sync is unchanched
  sys(4) = 1;     % -------------- '' --------------- , sync is (still) equal to 1
  %sys(5) = abs(-(4/u(3))*(t-round(t/u(3))*u(3)))-1;
  
elseif k == 10  %t_a_u_off    
  if status == 0  % -->Negative Saturation
      sys(1)=0;
      sys(2)=1;
  elseif status ==1 %-->Positive Saturation
      sys(1)=1;
      sys(2)=0;
  else
      sys(1)=0;
      sys(2)=0;
  end
  sys(3) = t_sync;% ---------- '' ----------  , t_sync is unchanched
  sys(4) = 1;     % ---------- '' ----------  , sync is (still) equal to 1
  %sys(5) = abs(-(4/u(3))*(t-round(t/u(3))*u(3)))-1;
  
elseif k == 11  %t_a_l_on-dt_plot    
  if status == 0  % -->Negative Saturation
      sys(1)=0;
      sys(2)=1;
  elseif status ==1 %-->Positive Saturation
      sys(1)=1;
      sys(2)=0;
  else
      sys(1)=0;
      sys(2)=0;
  end
  sys(3) = t_sync;% -------------- '' -------------  , t_sync is unchanched
  sys(4) = 1;     % -------------- '' -------------  , sync is (still) equal to 1
  %sys(5) = abs(-(4/u(3))*(t-round(t/u(3))*u(3)))-1;
  
elseif k == 12 %t_a_l_on
  if status==1 %-->Positive Saturation
      sys(1)=1;
      sys(2)=0;
  else
      sys(1)=0;
      sys(2)=1;      
  end
  sys(3) = t_sync;% ---------- '' ---------- , t_sync is unchanched
  sys(4) = 1;     % ---------- '' ---------- , sync is (still) equal to 1
  %sys(5) = abs(-(4/u(3))*(t-round(t/u(3))*u(3)))-1;
  
elseif k == 13 %t_next - dt_plot       
  if status==1 %-->Positive Saturation
      sys(1)=1;
      sys(2)=0;
  else
      sys(1)=0;
      sys(2)=1;      
  end
  sys(3) = t_sync;% ------------- '' ------------- , t_sync is unchanched
  sys(4) = 1;     % ------------- '' ------------- , sync is (still) equal to 1
  %sys(5) = abs(-(4/u(3))*(t-round(t/u(3))*u(3)))-1;
  
elseif k == 14;        
  sys(1) = 0;     % enable is turned off, da_u is set equal to 0   
  sys(2) = 0;     % -------- '' ------- , da_l is set equal to 0
  sys(3) = t_sync;% -------- '' ------- , t_sync is set on the output
  sys(4) = 0;     % -------- '' ------- , sync is set equal to 0
  %sys(5) = 0;
end

% end mdlOutputs

%
%=============================================================================
% mdlGetTimeOfNextVarHit
% Return the time of the next hit for this block.  Note that the result is
% absolute time.
%=============================================================================
%
function sys=mdlGetTimeOfNextVarHit(t,x,u,theta_phase,t_dead,t_min,type_signal,asym)

variable_vect = get_param(gcbh,'UserData');
dt_a_u_on = variable_vect(13);
U_st = u(1);
ksi_s=u(2);
T_tri = u(3);
enable = u(4);
k = variable_vect(1);
dt_plot = variable_vect(2);
t_sync = variable_vect(3);
t_next = variable_vect(4);
t_xtra = variable_vect(5);
t_a_u_on = variable_vect(6);
t_a_u_off = variable_vect(7);
t_a_l_on = variable_vect(8);
t_a_l_off = variable_vect(9);
u_st_a    = variable_vect(10);
status  = variable_vect(11);
status_  = variable_vect(12);
ignore_xtra  = variable_vect(14);
% triangle
%u_tri=abs(-(4/T_tri)*(t-round(t/T_tri)*T_tri))-1;
%u_tri=abs(-(4/T_tri)*(t_sync-round(t_sync/T_tri)*T_tri))-1;

if (k == 13 || k ==14 || enable==0)
  k = 0;    % re-initalization of k
else
  k = k+1;  % counter variable adds one
end % end if

if k == 0
  if type_signal == 2
     u_st_a =U_st;
  else    % type_signal == 1 or any other value 
     u_st_a =U_st;                 % "default" 
  end % end if
  
  if u_st_a > 1.0
    u_st_a =1.0;
  elseif u_st_a < -1.0    
    u_st_a =-1.0;   
  end
      
  dt_a_u_on = 0.25*(1-u_st_a)*T_tri;    
  status_= status;
  ignore_xtra=0;
  if dt_a_u_on < (t_min/2 + t_dead) && status == 1; %Positive to positive saturation
      dt_a_u_on=t_min/2+t_dead;
      dt_a_u_off=T_tri-t_min-t_dead;
      status = 1;
  elseif dt_a_u_on < (t_min+t_dead) && status == 1; %Positive to no saturation
      dt_a_u_off = T_tri-t_min/2;
      dt_a_u_on = t_min;
      status = 2;
      ignore_xtra = 1; %Ignore extra switching
  elseif dt_a_u_on < (T_tri-t_min-t_dead)/2 && status == 1; %Positive to no saturation
      dt_a_u_on = dt_a_u_on;
      dt_a_u_off = T_tri-dt_a_u_on;
      status = 2;
  elseif status == 1 %Positive to negative saturation
      dt_a_u_on = t_min*2;
      dt_a_u_off=T_tri-t_min; 
      status = 0;
  elseif dt_a_u_on < (t_min/2+t_dead) && status == 2; %No saturation to positive saturation
      dt_a_u_on = t_min/2+t_dead;
      dt_a_u_off = T_tri-dt_a_u_on;
      status=1;
  elseif dt_a_u_on < (T_tri - t_min - t_dead)/2 && status == 2; %No saturation to no saturation
      dt_a_u_on = dt_a_u_on;
      dt_a_u_off = T_tri - dt_a_u_on;
      status=2;
  elseif status == 2; %No saturation to negative saturation
      dt_a_u_on = t_min*2;
      dt_a_u_off = T_tri - t_min;
      status=0;
  elseif dt_a_u_on < t_min/2 && status == 0; %Negative to positive saturation;
      dt_a_u_on = t_min/2;
      dt_a_u_off = T_tri - t_min;
      status=1;
  elseif dt_a_u_on < (T_tri - t_min - t_dead)/2 && status == 0; %Negative to no saturation
      dt_a_u_on = dt_a_u_on;
      dt_a_u_off = T_tri - dt_a_u_on;
      status=2;
  elseif status == 0; %Negative to negative saturation
      dt_a_u_on = t_min;
      dt_a_u_off = T_tri - t_min;
      status=0;         
  end
     
  t_a_u_on = t+dt_a_u_on+t_dead;    % t_dead is the blanking time
  t_a_u_off = t+dt_a_u_off;
  t_a_l_on = t+dt_a_u_off+t_dead;    % t_dead is the blanking time
  t_a_l_off = t+dt_a_u_on;
  
  t_sync = t+T_tri/2;
  t_next = t+T_tri;
  t_xtra = t+t_dead;
  sys = t_xtra-dt_plot;
  
  if enable == 0    % checking if enable is off
    sys = t_next;
    k = 13; 
  end % end if  

elseif k == 1 %Xtra plot switching when going from higher band
  sys = t_xtra;
elseif k == 2 %Xtra switching when going from higher band
  sys = t_a_l_off-dt_plot;
elseif k == 3;
  sys = t_a_l_off;
elseif k == 4;
  sys = t_a_u_on-dt_plot;
elseif k ==5;
  sys = t_a_u_on;
elseif k == 6;
  sys = t_sync-dt_plot;
elseif k == 7;   
  sys = t_sync;
  
  if asym ==1
    if type_signal == 2
       u_st_a =U_st;
    else    % type_signal == 1 or any other value 
       u_st_a =U_st;                 % "default" 
    end % end if
    if u_st_a > 1.0
       u_st_a =1.0;
    elseif u_st_a < -1.0    
       u_st_a =-1.0;   
    end     
    dt_a_u_on = 0.25*(1-u_st_a)*T_tri;    
    status_= status;
    ignore_xtra=0;
    if dt_a_u_on < (t_min/2 + t_dead) && status == 1; %Positive to positive saturation
        dt_a_u_off=T_tri-t_min-t_dead;
        status = 1;
    elseif dt_a_u_on < (t_min+t_dead) && status == 1; %Positive to no saturation
        dt_a_u_off = T_tri-t_min/2;
        status = 2;
        ignore_xtra = 1; %Ignore extra switching
    elseif dt_a_u_on < (T_tri-t_min-t_dead)/2 && status == 1; %Positive to no saturation
        dt_a_u_off = T_tri-dt_a_u_on;
        status = 2;
    elseif status == 1 %Positive to negative saturation
        dt_a_u_off=T_tri-t_min; 
        status = 0;
    elseif dt_a_u_on < (t_min/2+t_dead) && status == 2; %No saturation to positive saturation
        dt_a_u_on = t_min/2+t_dead;
        dt_a_u_off = T_tri-dt_a_u_on;
        status=1;
    elseif dt_a_u_on < (T_tri - t_min - t_dead)/2 && status == 2; %No saturation to no saturation
        dt_a_u_off = T_tri - dt_a_u_on;
        status=2;
    elseif status == 2; %No saturation to negative saturation
        dt_a_u_off = T_tri - t_min;
        status=0;
    elseif dt_a_u_on < t_min/2 && status == 0; %Negative to positive saturation;
        dt_a_u_off = T_tri - t_min;
        status=1;
    elseif dt_a_u_on < (T_tri - t_min - t_dead)/2 && status == 0; %Negative to no saturation
        dt_a_u_off = T_tri - dt_a_u_on;
        status=2;
    elseif status == 0; %Negative to negative saturation
        dt_a_u_off = T_tri - t_min;
        status=0;         
    end    
    t_a_u_off = t-T_tri/2+dt_a_u_off;
    t_a_l_on = t-T_tri/2+dt_a_u_off+t_dead;    % t_dead is the blanking time
  end % End asym
            
elseif k == 8;  
  sys = t_a_u_off-dt_plot; 
elseif k == 9;  
  sys = t_a_u_off; 
elseif k == 10;  
  sys = t_a_l_on-dt_plot;
elseif k == 11;  
  sys = t_a_l_on;
elseif k == 12;  
  sys = t_next-dt_plot; 
elseif k == 13;  
  sys = t_next;
else    
  k = -1;    % re-initalization of k 
end    % end if
set_param(gcbh,'UserData',[k,dt_plot,t_sync,t_next,t_xtra,t_a_u_on,t_a_u_off,t_a_l_on,t_a_l_off,u_st_a,status,status_,dt_a_u_on,ignore_xtra]);
% writes the variable values to the block's workspace memory

% end mdlGetTimeOfNextVarHit

%
%=============================================================================
% mdlTerminate
% Perform any end of simulation tasks.
%=============================================================================
%
function sys=mdlTerminate(t,x,u)

sys = [];

% end mdlTerminate
