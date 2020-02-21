function [preamble, header, psdu] = generatePacket(data, dataRate)
    % this packet frame generation for PHY layer is only for 
    % 'long' PLCP preamble.
    % 
    % * PREAMBLE and HEADER fields shall be transmitted using the 1 Mbps DBPSK
    % modulation (using Barker Sequence)
    % * SIGNAL and SERVICE fields indicate the rate and the modulation type,
    % respectively.
    
    % type = 'long';
    if nargin < 2
        error('[Error] Not enough input parameters!')
    end
    
    % initialization (row vectors)
    psdu = data;
    
    % -------PREAMBLE-------
    % generate scrambled SYNC bits
    orig = ones([128,1]);
    scramInit = [1 1 0 1 1 0 0]';
    syncScrambled = wlanScramble(orig,scramInit); % MSB to LSB
    
    % generate SFD bits
    sfd = fliplr([1 1 1 1 0 0 1 1 1 0 1 0 0 0 0 0]); % LSB to MSB
    preamble = [syncScrambled', sfd]; % row vector
    
    % -------HEADER-------
    % generate SIGNAL
    signal = []; % MSB to LSB
    switch dataRate
        case 1
            signal = [zeros(1,4), 1 0 1 0]; % X'0A'
        case 2
            signal = [zeros(1,3), 1 0 1 0 0]; % X'14'
        case 5.5
            signal = [zeros(1,2), 1 1 0 1 1 1]; % X'37'
        case 11
            signal = [zeros(1,1), 1 1 0 1 1 1 0]; % X'6E'
    end
    
    % generate SERVICE bits
    reserved = 0;
    clocksbit = 1; % locked
    modulation = 0; % cck = 0, pbcc = 1
    lengthExt = 0;
    
    % for CCK -- 5.5 and 11 Mbps
    numTxOctets = size(reshape(data, 8, []),1);
    P = 0; % 0 for CCK
    temp = (numTxOctets + P) * 8 / dataRate; 
    tempCeiling = ceil(temp);
    % reference -- time: 8 microseconds per octet in PPDU
    length = de2bi(round(temp), 16, 'left-msb');
    if dataRate == 11 && (tempCeiling - temp)>=8/11
        lengthExt = 1;
    else
        lengthExt = 0;
    end
    % MSB to LSB
    service = [ reserved, reserved, clocksbit, modulation, ...
                reserved, reserved, reserved, lengthExt];
    
    % at the receiver side, number of octets in MPDU is calculated as:
    % nOctets = floor(((length * dataRate)/8)-P)-lengthExt;
    
    % generate CCITT CRC-16 frame check sequence (FCS) / MSB to LSB
    input = [signal, service, lengthExt];
    % 1. preset to all one's
    % 2. shift signal,service,length fields through the shift register
    % 3. take one's complement of the remainder
    % 4. transmit out serial X^15 first
    % construct a CRC generator with a polynomial defined by
    
    % 1 + x^5 + x^12 + x^16
    poly = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 1];
    init = ones([1,15]);
    h = crc.generator('Polynomial',poly, ...
                        'InitialState',init);
    crcEncoded = generate(h, input');
    crcOnesComp = double(not(num2str(crcEncoded) - '0'));  % one's complement
    crcBits = fliplr(crcOnesComp(end-15:end)');
    
    % put together header
    header = [signal, service, length, crcBits];
end    
    
    