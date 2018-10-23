%%%%%%%%%%%%%%%%%%%%%%% PREMABLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Author:          Tumisang Kalagobe
%% Date completed:  22/08/2018
%% Course:          Mechanical Engineering Design (MECN4005)
%% Institution:     University of the Witwatersrand
%% Description:     Design of an energy recovery system for cryptocurrency
%%                  mining activities. This program is dedicated to
%%                  evaluating the condenser of a small organic Rankine 
%%                  cycle.
%%                  The heat transfer coefficient is evaluated to within
%%                  30% of the assumed value.

clc
clear

tic 
%%%%%%%%%%%%%%%%%%%%%%% MAIN CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 10; % number of points of interest
min_W_net = 0.065; % 5% net work that must be recovered (kW)
max_W_net = 0.1345; % 10% net work that must be recovered (kW)
W_net = linspace(min_W_net, max_W_net, n);
eta = 0.1; % thermal efficiency (W_net/Q_in)
Q_in = W_net./eta; % heat energy input (kW)
cp_r = 1.36 ; % refrigerant specific heat (kJ/kg.K)
cp_a = 1.003; % air specific heat (kJ/kg.K)
L = 0.3; % length of a single tube (m)
di = 0.005; % inner diameter of the pipe (m)
t = di*0.2; % combined thickness of tubes - i.e t*2 (m)
do = di + t; % external pipe diameter (m)

U = 0; % assumed heat transfer coefficient (kW/m^2.K)
tol = 1; % percentage tolerance 
z = 0;
maximum_error = 30;
while maximum_error > tol
    
    [n_tubes, m_dot_a, m_dot_r] = effNTU(n, Q_in, cp_a, cp_r, L, do, U);
    [error, U_calc] = coeffError(m_dot_r, m_dot_a, do, di, n, U, L, cp_r);
    error = 100*error;
    maximum_error = max(abs(error));
    required_tubes = max(n_tubes);
    U = U + 0.001;
    z = z + 1;
    
    if U > 10
        tol = tol + 1;
        U = 0.05;
    elseif tol > 30
        break
    end
end 

if maximum_error <= 30
    disp("Error is within 30%!")
    disp(" ")
    
    error_statement = "Error in HT coefficient = " + string(maximum_error);
    disp(error_statement)
    disp(" ")
    
    tube_statement = "No. of tubes required = " + string(required_tubes);
    disp(tube_statement)
    disp(" ")
else 
    disp("Error is greater than 30%")
    disp("Please check parameters.")
end

toc
%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n_tubes, m_dot_a, m_dot_r] = effNTU(n, Q_in, cp_a, cp_r, L, do, U)

%% ------------------------- Description ----------------------------------
%%  Estimating the required mass flow rates of air and R245fa as well as 
%%  the number of tubes required to meet the energy requirements via the
%%  effectiveness-NTU method for a single tube, single pass 
%%  shell and tube HX.

%% Variables:
%%  T_ref_out    =  refrigerent outlet temperature (K)
%%  T_ref_in     =  refrigerent inlet temperature (K)
%%  T_air_in     =  air inlet temperature (K)
%%  T_air_out    =  air outlet temperature (K)
%%  Q_max        =  maximum possible heat transfer in the system (kW)
%%  epsilon      =  effectiveness 
%%  NTU          =  Net Transfer Units
%%  A_evap_tubes =  Required tube area (m^2)

%% Outputs:
%%  n_tubes      =  Number of tubes required
%%  m_dot_r      =  refrigerant mass flow rate (kg/s)
%%  m_dot_a      =  air mass flow rate (kg/s)
%% ------------------------------------------------------------------------
    T_air_in = 353;
    T_ref_out = 340;
    T_ref_in = 303; 
    T_air_out = T_ref_out + 5;
    dTr = T_ref_out - T_ref_in;
    dTa = T_air_in - T_air_out;
    m_dot_r = Q_in./(cp_r.*dTr); 
    m_dot_a = Q_in./(cp_a.*dTa);

    % effectiveness-NTU approach
    C_air = cp_a.*m_dot_a;
    C_ref = cp_r.*m_dot_r.*ones(1,n);
    C_min = zeros(1,n);
    C_max = zeros(1,n);
    for i=1:1:n
        if C_air(i) > C_ref(i)
            C_min(i) = C_ref(i);
            C_max(i) = C_air(i);
        else 
            C_min(i) = C_air(i);
            C_max(i) = C_ref(i);
        end
    end

    Cr = C_min./C_max;
    Q_max = C_min.*(T_air_in - T_ref_in);
    epsilon = Q_in./Q_max;
    
    NTU = zeros(1,n);
    E = zeros(1,n);
    for i=1:1:n
        E(i) = (2./epsilon(i) - 1 + Cr(i))./(sqrt(1 + Cr(i)^2));
        NTU(i) = -sqrt(1 + Cr(i)^2)*log((E(i) - 1)./(E(i) + 1));
    end

    % Geometry evaluation
    A_evap_tubes = NTU.*C_min./U;
    n_tubes = round(A_evap_tubes./(pi.*do.*L),0);
end

function [error, U_calc] = coeffError(m_dot_r, m_dot_a, do, di, n, U, L, cp_r)
%% ------------------------- Description ----------------------------------
%%  Determining the overall heat transfer coefficient of the HX by use of 
%%  standard HT relations for a tube in crossflow. The Sieder-Tate HT 
%%  correlation is used to evaluate the internal tube convective HT 
%%  coefficient.
%%
%% Variables:
%%  A_inlet_air     =   air inlet area (m^2)
%%  miu             =   dynamic viscosity (kg/m.s)
%%  Pr              =   Prandtl number 
%%  h               =   convective heat transfer coefficient (kW/m^2.K)
%%  k               =   conductive heat transfer coefficient (kW/m^2.K)
%%
%% Outputs:
%%  error           =   percentage error between the assumed and calculated
%%                      heat transfer coeficcient values
%%  U_calc          =   calculated value for the overall heat transfer
%%                      coefficient
%% ------------------------------------------------------------------------
    %% external heat transfer - cylinder in crossflow
    A_inlet_air = (pi*(0.2*L)^2)/4;
    miu_air = 208e-7;
    k_air = 30e-6; 
    Pr_air = 0.7;
    Re_air = m_dot_a.*do./(A_inlet_air*miu_air);
    h_ext = zeros(1,n);
    for i = 1:1:n
        if Re_air(i) < 40
            C = 0.911;
            m = 0.385;
        elseif Re_air(i) >= 40 && Re_air(i) < 4000
            C = 0.683;
            m = 0.466;
        elseif Re_air(i) >= 4000 && Re_air(i) < 40000
            C = 0.193;
            m = 0.618;
        elseif Re_air(i) >= 40000 && Re_air(i) < 400000
            C = 0.027;
            m = 0.805;
        end
        h_ext(i) = (k_air/do)*(C.*(Re_air(i)^m)*(Pr_air^0.333));
    end
    
    %% internal flow - Sieder-Tate correlation
    k_tube = 15e-3;
    miu_ref = 10.69e-6;
    miu_ref_w = 17.58e-6;
    k_ref = 13.87e-6;
    Pr_ref = cp_r*miu_ref/k_ref;
    A_tube_cross_section = pi*(di^2)/4;
    Re_ref = (m_dot_r.*di)./(A_tube_cross_section.*miu_ref);
    PrMiu = (Pr_ref^(1/3)).*(miu_ref/miu_ref_w)^0.14;
    h_int = (k_ref/di).*(0.023.*(Re_ref.^0.8).*PrMiu);

    %% Overall HT coefficient
    internal_conv = 1./h_int;
    external_conv = (di)./(do*h_ext);
    tube_cond = (di*log(do/di))/(2*k_tube);
    U_calc = (internal_conv + tube_cond + external_conv).^(-1); % W/m^2.K
    error = (U_calc - U)./U; % must be between 0 and 30%
end 
