classdef FoilParamsESC < handle
    properties
        esc_type string % GB / NB
        f               % dither frequency
        A               % dither magnitude
        phase_A         % dither phase
        fhp             % cutoff frequency HPF
        % flp             % cutoff frequency LPF
        K               % integrator gain (>0 maximum)
        B               % demodulation magnitude
        phase_B         % demodulation phase
        % NB-specific
        wric            % ricatti filter param
        ric_0           % ricatti filter initial estimate
    end
    
    methods
        function obj = Foil(esc_type, f, A, phase_A, fhp, K, B, phase_B, wric, ric_0) % Constructor method
            if strcmp(esc_type, "NBESC") && nargin < 9
                disp('Error: Not enough arguments for NBESC class.\n')
            end
            obj.f       = f;
            obj.A       = A;
            obj.phase_A = phase_A;
            obj.fhp     = fhp;
            obj.K       = K;
            obj.B       = B;
            obj.phase_B = phase_B;

            if strcmp(esc_type, "NBESC") 
                obj.wric  = wric;
                obj.ric_0 = ric_0;
            end
        end
    end
end

    