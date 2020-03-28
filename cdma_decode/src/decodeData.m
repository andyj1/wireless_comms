function [result] = decodeData(data)
    
    % takes data signal only (w/o pilot symbols) and demodulates 
    
    % make column vector
    if size(data,2) > 1
        data = data';
    end
    
    % decodes BPSK modulated signal into -1, 0, 1
    % such that fo reach bit,
    % -1 if x < -0.5, 0 if -0.5 <= x <= 0.5, and 1 if x > 0.5
    result = zeros(size(data));
    for i = 1:length(data)
        if data(i) > 0.5
            result(i) = 1;
        elseif data(i) >= -0.5 && data(i) <= 0.5
            result(i) = 0;
        elseif data(i) < -0.5
            result(i) = -1;
        end
    end

end