function [result] = decodePilot(pilot_hadamard, PNsequence)
    
    % decodes pilot using Walsh channel 0

    % Walsh code
    Hadamard = pilot_hadamard(1:length(PNsequence));
    % BPSK demodulate
    Hadamard_demod = bpsk.demodBPSK(Hadamard);
    % despread
    result = xor(Hadamard_demod, PNsequence);

end