classdef ManipulateData
    % Class that manage multiple functions


    properties(Constant)
        
    end
    
    
    methods(Static = true)
        function Ke = ComputeKrotation(E, Ic, H, K_factor)
            % Rotational stiffness
            if nargin < 4
                K_factor = 3; % Cantelever
            end
            Ke = K_factor*E*Ic/H;
        end

        function Ks = ComputeKshear(d, bf, tf, tw, A, H, G, K_shear_factor)
            % Rotational stiffness
            if nargin < 8
                K_shear_factor = 1; % Cantelever
            end
            Ks = G*A*(K_shear_factor*H) / (0.85+2.32*bf*tf/d/tw); % (Charney et al. 2005) 
        end

        function Kc = ComputeKCol(d, bf, tf, tw, A, Ic, H, E, H_B, H_T)
            % H_B = length column Bottom without panel zone
            % H_T = length column Top without panel zone

            G = E/2.6;
            % Bottom
            Kb_B = ManipulateData.ComputeKrotation(E, Ic, H_B);
            Ks_B = ManipulateData.ComputeKshear(d, bf, tf, tw, A, H_B, G);
            KrotB = (Kb_B*Ks_B)/(Kb_B+Ks_B); % kNmm - Column rotational stiffness (cantilevel moment at base - chord rotation)
            KeB = KrotB/H_B/H_B;

            % Top
            Kb_T = ManipulateData.ComputeKrotation(E, Ic, H_T);
            Ks_T = ManipulateData.ComputeKshear(d, bf, tf, tw, A, H_T, G);
            KrotT = (Kb_T*Ks_T)/(Kb_T+Ks_T); % kNmm - Column rotational stiffness (cantilevel moment at base - chord rotation)
            KeT = KrotT/H_T/H_T;

            Kc = (1/KeB + 1/KeT)^(-1)*H;
        end
        
        function theta_adjusted = AdjustThetaWithn(My_star, K, theta, n)
            % Adjust the stiffness of the elastic regime with n and
            % translate the rest accordingly
            
            if nargin < 4
                n = 10;
            end
            
            theta_y = My_star/K;
            theta_adjusted = theta;

            % Use n
            i = 1;
            while i < length(theta_adjusted) && theta_adjusted(i) <= theta_y
                theta_adjusted(i) = theta_adjusted(i)*(n+1);
                i = i + 1;
            end
            theta_adjusted(i:end) = theta_adjusted(i:end)+theta_y*(n+1);
        end
    end
    
    
end
