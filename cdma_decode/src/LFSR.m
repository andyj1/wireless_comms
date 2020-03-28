function [PNsequence] = LFSR(poly, initial_state)

    % feedback taps of order m
    % An LFSR of any given size m (number of register)
    % is capable of producing every possible state during 
    % the period N = 2^m - 1 shifts
    % --> maximal length sequance (m-sequence)
    m = max(poly); % highest order
    
    % initialize PNsequence (binary) to the initial state
    states = zeros([m, 1]);
    states(initial_state) = 1;

    % m-bit register produces an m-sequence of period N (=2^m - 1)
    N = 2^m - 1;   % total sequence length
    PNsequence = zeros([N 1]);
    for i = 1:N
        PNsequence(i) = states(end);
        temp = mod(sum(states(poly)),2);
        states(2:end) = states(1:end-1);
        states(1) = temp;
    end
    
end
