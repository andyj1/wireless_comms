classdef bpsk
    
    methods(Static)
        function [result] = modBPSK(data)
            
            % modulates binary bits to BPSK 
            % such that 1 -> 1, 0 -> -1
            
            % make column vector
            if size(data,2) > 1
                data = data';
            end


            % initialize modulator
            bpskModulator = comm.BPSKModulator;
            bpskModulator.PhaseOffset = 0;

            % make column vector
            if size(data,2) > 1
                data = data';
            end

            % BPSK modulation, but reverse
            result = real(-1*bpskModulator(data));

        end 
        
        function [result] = demodBPSK(data)

            % demodulates BPSK to binary bits
            % such that 1 -> 1, -1 -> 0
            
            % make column vector
            if size(data,2) > 1
                data = data';
            end

            result = double(data > 0);

        end
        
        
    end
    
end