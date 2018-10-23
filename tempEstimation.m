%%%%%%%%%%%%%%%%%%%%%%% PREMABLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Author:          Tumisang Kalagobe
%% Date completed:  22/08/2018
%% Description:     First pass estimation of the required temperatures at
%%                  various inlet conditions for both the refrigerant and 
%%                  air at different heat inputs.

clc
clear

%%%%%%%%%%%%%%%%%%%%%%% MAIN FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 4; % number of points of interest
eta = 0.1; % thermal efficiency (W_net/Q_in)
W_net = linspace(0.0675,0.1375,n); % net work recovered in the system (kW)
Q_in = W_net./eta; % heat energy input (kW)
cp_r = 1.36 ; % refrigerant specific heat (kJ/kg.K)
cp_a = 1.003; % air specific heat (kJ/kg.K)
L = 0.3; % length of a single tube (m)
d_tube = 0.005; % diameter of the pipe (m)
U = 0.5; % assumed heat transfer coefficient (kW/m^2.K)

main(n, Q_in, cp_a, cp_r, L, d_tube, U)

%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main(n, Q_in, cp_a, cp_r, L, d_tube, U)

% comparing the effects of varying refrigerant output temperatures
compare_T_ref_out()

% comparing the effects of varying air input temperatures
compare_T_air_in()

% comparing the effects of varying refrigerant input temperatures
compare_T_ref_in()

    function compare_T_ref_out()
        % crossflow run
        for run = 1:1:n
            T_air_in = 353;
            T_ref_out = linspace(310,343,n);
            T_ref_in = 303; 
            T_air_out = T_ref_out(run) + 5;
            dTr = T_ref_out(run) - T_ref_in;
            dTa = T_air_in - T_air_out;
            m_dot_r = Q_in./(cp_r.*dTr); 
            m_dot_a = Q_in./(cp_a.*dTa);

            % effectiveness-NTU approach
            C_a = cp_a.*m_dot_a;
            C_R = cp_r.*m_dot_r;
            C_min = zeros(1,n);
            C_max = zeros(1,n);
            for i=1:1:length(C_a)
                if C_a(i) > C_R(i)
                    C_min(i) = C_R(i);
                    C_max(i) = C_a(i);
                else 
                    C_min(i) = C_a(i);
                    C_max(i) = C_R(i);
                end
            end

            Cr = C_min./C_max;
            Q_max = C_min.*(T_air_in - T_ref_in);
            epsilon = Q_in./Q_max; % effectiveness
            NTU = zeros(1,n);

            for i=1:1:length(Cr) %% crossflow single pass
                if C_min == C_a(i)
                    NTU(i) = -(1./Cr(i))*log(Cr(i).*log(1-epsilon(i)) + 1);
                else 
                    NTU(i) = -log(1 + (1./Cr(i)).*log(1 - epsilon(i).*Cr(i)));
                end
            end

            % heat transfer coefficient wrt NTU and area
            A_evap_tubes = NTU.*C_min./U;
            n_tubes = round(A_evap_tubes./(pi.*d_tube.*L),0);

            % plots
            legend_string = "T_{R245fa,out,CF} = " + string(round(T_ref_out(run),1)) + " K";
            if run == 1
                point = 'o';
            elseif run == 2
                point = '*';
            elseif run == 3
                point = '+';
            elseif run == 4
                point = '<';
            end

            figure(1)
            plot(Q_in.*1000,n_tubes, point)
            xlabel('Q_{in} (W)')
            ylabel('Number of tubes required')
            legend(legend_string)
            legend('boxoff')
            legend('Location','northeastoutside')
            legend('AutoUpdate','off')
            title('Varying R245fa outlet temperature')
            hold on
        end %% end of crossflow 

        % Shell and tube run 
        for run = 1:1:n
            T_air_in = 353;
            T_ref_out = linspace(310,343,n);
            T_ref_in = 303; 
            T_air_out = T_ref_out(run) + 5;
            dTr = T_ref_out(run) - T_ref_in;
            dTa = T_air_in - T_air_out;
            m_dot_r = Q_in./(cp_r.*dTr); 
            m_dot_a = Q_in./(cp_a.*dTa);

            % effectiveness-NTU approach
            C_a = cp_a.*m_dot_a;
            C_R = cp_r.*m_dot_r;
            C_min = zeros(1,n);
            C_max = zeros(1,n);
            for i=1:1:length(C_a)
                if C_a(i) > C_R(i)
                    C_min(i) = C_R(i);
                    C_max(i) = C_a(i);
                else 
                    C_min(i) = C_a(i);
                    C_max(i) = C_R(i);
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

            % heat transfer coefficient wrt NTU and area
            A_evap_tubes = NTU.*C_min./U;
            n_tubes = round(A_evap_tubes./(pi.*d_tube.*L),0);

            % plots
            output_temp = T_ref_out;
            legend_string = "T_{R245fa,out,ST} = " + string(round(output_temp,1)) + " K";
            if run == 1
                point = '-';
            elseif run == 2
                point = '--';
            elseif run == 3
                point = ':';
            elseif run == 4
                point = '-.';
            end

            figure(1)
            plot(Q_in.*1000,n_tubes, point)
            xlabel('Q_{in} (W)')
            ylabel('Number of tubes required')
            legend(legend_string)
            legend('boxoff')
            legend('Location','northeastoutside')
            hold on
        end % end of shell and tube
        %saveas(gcf,'T_R245fa_out.pdf')
    end 

    function compare_T_air_in()
        % crossflow run
        for run = 1:1:n        
            T_air_in = linspace(345,353,n);
            T_ref_out = 340;
            T_ref_in = 303; 
            T_air_out = T_ref_out + 5;
            dTr = T_ref_out - T_ref_in;
            dTa = T_air_in(run) - T_air_out;
            m_dot_r = Q_in./(cp_r.*dTr); 
            m_dot_a = Q_in./(cp_a.*dTa);

            % effectiveness-NTU approach
            C_a = cp_a.*m_dot_a;
            C_R = cp_r.*m_dot_r;
            C_min = zeros(1,n);
            C_max = zeros(1,n);
            for i=1:1:n
                if C_a(i) > C_R(i)
                    C_min(i) = C_R(i);
                    C_max(i) = C_a(i);
                else 
                    C_min(i) = C_a(i);
                    C_max(i) = C_R(i);
                end
            end

            Cr = C_min./C_max;
            Q_max = C_min.*(T_air_in(run) - T_ref_in);
            epsilon = Q_in./Q_max; % effectiveness
            NTU = zeros(1,n);
            for i=1:1:length(Cr) %% crossflow single pass
                if C_min == C_a(i)
                    NTU(i) = -(1./Cr(i))*log(Cr(i).*log(1-epsilon(i)) + 1);
                else 
                    NTU(i) = -log(1 + (1./Cr(i)).*log(1 - epsilon(i).*Cr(i)));
                end
            end

            % heat transfer coefficient wrt NTU and area
            A_evap_tubes = NTU.*C_min./U;
            n_tubes = round(A_evap_tubes./(pi.*d_tube.*L),0);

            % plots
            legend_string = "T_{air,in,CF} = " + string(round(T_air_in,1)) + " K";
            if run == 1
                point = 'o';
            elseif run == 2
                point = '*';
            elseif run == 3
                point = '+';
            elseif run == 4
                point = '<';
            elseif run == 5
                point = 'x';
            end

            figure(3)
            plot(Q_in.*1000,n_tubes, point)
            xlabel('Q_{in} (W)')
            ylabel('Number of tubes required')
            legend(legend_string)
            legend('boxoff')
            legend('Location','northeastoutside')
            title('Varying air inlet temperatures')
            hold on
        end 

        % shell and tube run
        for run = 1:1:n
            T_air_in = linspace(345,353,n);
            T_ref_out = 340;
            T_ref_in = 303; 
            T_air_out = T_ref_out + 5;
            dTr = T_ref_out - T_ref_in;
            dTa = T_air_in(run) - T_air_out;
            m_dot_r = Q_in./(cp_r.*dTr); 
            m_dot_a = Q_in./(cp_a.*dTa);

            % effectiveness-NTU approach
            C_a = cp_a.*m_dot_a;
            C_R = cp_r.*m_dot_r;
            C_min = zeros(1,n);
            C_max = zeros(1,n);
            for i=1:1:n
                if C_a(i) > C_R(i)
                    C_min(i) = C_R(i);
                    C_max(i) = C_a(i);
                else 
                    C_min(i) = C_a(i);
                    C_max(i) = C_R(i);
                end
            end

            Cr = C_min./C_max;
            Q_max = C_min.*(T_air_in(run) - T_ref_in);
            epsilon = Q_in./Q_max; % effectiveness
            NTU = zeros(1,n);
            E = zeros(1,n);
            for i=1:1:length(Cr) %% shell and tube
                E(i) = (2./epsilon(i) - 1 + Cr(i))./(sqrt(1 + Cr(i)^2));
                NTU(i) = -sqrt(1 + Cr(i)^2)*log((E(i) - 1)./(E(i) + 1));
            end

            % heat transfer coefficient wrt NTU and area
            A_evap_tubes = NTU.*C_min./U;
            n_tubes = round(A_evap_tubes./(pi.*d_tube.*L));

            % plots
            legend_string = "T_{air,in,ST} = " + string(round(T_air_in,1)) + " K";
            if run == 1
                point = '-';
            elseif run == 2
                point = '--';
            elseif run == 3
                point = ':';
            elseif run == 4
                point = '-.';
            end

            figure(3)
            plot(Q_in.*1000,n_tubes, point)
            xlabel('Q_{in} (W)')
            ylabel('Number of tubes required')
            legend(legend_string)
            legend('boxoff')
            legend('Location','northeastoutside')
            hold on
        end 
        %saveas(gcf,'T_air_in (crossflow).pdf')
    end

    function compare_T_ref_in()
        % crossflow run
        for run = 1:1:n
            T_air_in = 353;
            T_ref_out = 340;
            T_ref_in = linspace(298,320,n); 
            T_air_out = T_ref_out + 5;
            dTr = T_ref_out - T_ref_in(run);
            dTa = T_air_in - T_air_out;
            m_dot_r = Q_in./(cp_r.*dTr); 
            m_dot_a = Q_in./(cp_a.*dTa);

            % effectiveness-NTU approach
            C_a = cp_a.*m_dot_a;
            C_R = cp_r.*m_dot_r;
            C_min = zeros(1,n);
            C_max = zeros(1,n);
            for i=1:1:n
                if C_a(i) > C_R(i)
                    C_min(i) = C_R(i);
                    C_max(i) = C_a(i);
                else 
                    C_min(i) = C_a(i);
                    C_max(i) = C_R(i);
                end
            end

            Cr = C_min./C_max;
            Q_max = C_min.*(T_air_in - T_ref_in(run));
            epsilon = Q_in./Q_max; % effectiveness
            NTU = zeros(1,n);
            for i=1:1:n %% crossflow single pass
                if C_min == C_a(i)
                    NTU(i) = -(1./Cr(i))*log(Cr(i).*log(1-epsilon(i)) + 1);
                else 
                    NTU(i) = -log(1 + (1./Cr(i)).*log(1 - epsilon(i).*Cr(i)));
                end
            end

            % heat transfer coefficient wrt NTU and area
            A_evap_tubes = NTU.*C_min./U;
            n_tubes = round(A_evap_tubes./(pi.*d_tube.*L),0);

            % plots
            legend_string = "T_{R245fa,in, CF} = " + string(round(T_ref_in,1)) + " K";
            if run == 1
                point = 'o';
            elseif run == 2
                point = '*';
            elseif run == 3
                point = '+';
            elseif run == 4
                point = '<';
            elseif run == 5
                point = 'x';
            end

            figure(5)
            plot(Q_in.*1000,n_tubes, point)
            xlabel('Q_{in} (W)')
            ylabel('Number of tubes required')
            legend(legend_string)
            legend('boxoff')
            legend('Location','northeastoutside')
            title('Number of tubes required in relation to heat input for a crossflow HX')
            hold on
        end 

        % shell and tube run
        for run = 1:1:n
            T_air_in = 353;
            T_ref_out = 340;
            T_ref_in = linspace(298,320,n); 
            T_air_out = T_ref_out + 5;
            dTr = T_ref_out - T_ref_in(run);
            dTa = T_air_in - T_air_out;
            m_dot_r = Q_in./(cp_r.*dTr); 
            m_dot_a = Q_in./(cp_a.*dTa);

            % effectiveness-NTU approach
            C_a = cp_a.*m_dot_a;
            C_R = cp_r.*m_dot_r;
            C_min = zeros(1,n);
            C_max = zeros(1,n);
            for i=1:1:n
                if C_a(i) > C_R(i)
                    C_min(i) = C_R(i);
                    C_max(i) = C_a(i);
                else 
                    C_min(i) = C_a(i);
                    C_max(i) = C_R(i);
                end
            end

            Cr = C_min./C_max;
            Q_max = C_min.*(T_air_in - T_ref_in(run));
            epsilon = Q_in./Q_max; % effectiveness
            NTU = zeros(1,n);
            E = zeros(1,n);
            for i=1:1:length(Cr) %% shell and tube
                E(i) = (2./epsilon(i) - 1 + Cr(i))./(sqrt(1 + Cr(i)^2));
                NTU(i) = -sqrt(1 + Cr(i)^2)*log((E(i) - 1)./(E(i) + 1));
            end

            % heat transfer coefficient wrt NTU and area
            A_evap_tubes = NTU.*C_min./U;
            n_tubes = round(A_evap_tubes./(pi.*d_tube.*L));

            % plots
            legend_string = "T_{R245fa,in,ST} = " + string(round(T_ref_in,1)) + " K";
            if run == 1
                point = '-';
            elseif run == 2
                point = '--';
            elseif run == 3
                point = ':';
            elseif run == 4
                point = '-.';
            end

            figure(5)
            plot(Q_in.*1000,n_tubes, point)
            xlabel('Q_{in} (W)')
            ylabel('Number of tubes required')
            legend(legend_string)
            legend('boxoff')
            legend('Location','northeastoutside')
            title('Varying refrigerant inlet temperatures')
            hold on
        end 
    end 

end
