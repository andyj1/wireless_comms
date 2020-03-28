function [result] = translate(data, PNsequence, Hchan)

    % converts decoded BPSK signal to ASCII code

    % 1. divide data into frame chunks
    % - frames1 : frames from 2nd to 2nd to last
    % - frames2 : last frame, without zeros
    frames = reshape(data, length(PNsequence), []);
    frames2 = frames(:, end); % in case last frame has fewer than 3 characters
    frames1 = frames(:, 1:end-1);
    
    % 2. despread using the PN sequence
    % - data1 : take 3 character-long chips (cpf: 64, # char: 3)
    cpf = 64; num_char = 3;
    data1 = zeros(cpf*num_char, size(frames1,2));
    for i = 1:size(frames1,2)
        frames1_demod = bpsk.demodBPSK( frames1(frames1(:,i)~=0,i) );
        % check sizes for xor operation
        % [size(frames1_demod), size(PNsequence(1:length(frames1_demod')))]
        data1(:,i) = xor( frames1_demod, PNsequence(1:cpf*num_char) );
    end
    % - data2 : only take nonzero symbols
    frames2_demod = bpsk.demodBPSK(frames2(frames2~=0));
    data2 = xor(frames2_demod, PNsequence(1:length(frames2_demod)));
    % concatenate data1 and data2
    data = [data1(:); data2(:)];
    
    % Walsh decode
    WalshDecoded = Hchan * reshape(bpsk.modBPSK(data),8,[]) / 8;
    % BPSK demodulate
    decoded = bpsk.demodBPSK(WalshDecoded);
    % convert to ASCII characters
    ascii_dec = bi2de(reshape(decoded,8,[]).');
    result = char(ascii_dec.');
    
end