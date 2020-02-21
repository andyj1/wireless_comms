function [checked] = checkFrame(preambleTx, headerTx, preambleRx, headerRx, dataRate)
    % checks the values for contents in the preamble and header
    if nargin < 4
        error('[Error] Not enough input parameters!')
    end
    
    valid = 1;
    checked = true;
    
    % -------PREAMBLE-------
    % check SYNC by descrambling to see if it returns 'orig'
    orig = ones([128,1]);
    scramInit = [1 1 0 1 1 0 0]';
    syncDescrambled = wlanScramble(preambleRx(1:128),scramInit);
    if syncDescrambled ~= orig
        error('sync not matching!')
%     else
%         disp('SYNC checked')
    end
    
    % check SFD to see if it returns 'orig'
    orig = fliplr([1 1 1 1 0 0 1 1 1 0 1 0 0 0 0 0]);
    sfd = preambleRx(129:129+16-1);
    if sfd ~= orig
        error('sfd not matching!')
%     else
%         disp('SFD checked')
    end
    
    % -------HEADER-------
    % check SIGNAL to see if it returns 'orig'
    signal = headerRx(1:8);
    orig = [];
    switch dataRate
        case 1
            orig = [zeros(1,4), 1 0 1 0];
        case 2
            orig = [zeros(1,3), 1 0 1 0 0];
        case 5.5
            orig = [zeros(1,2), 1 1 0 1 1 1];
        case 11
            orig = [zeros(1,1), 1 1 0 1 1 1 0];
    end
    if ~isequal(signal, orig)
        error('signal not matching!')
%     else
%         disp('signal checked')
    end
    
    % check SERVICE to see if it returns 'orig'
    service = headerRx(9:9+8-1);
    orig = 0; % CCK bit
    % check modulation scheme (CCK)
    if service(4) ~= orig   
        error('MCS is not CCK!')
%     else
%         disp('service bit checked');
    end
    
    % acquire LENGTHEXT bit
    lengthext = headerRx(17:17+16-1);
    
    % check 16-bit CRC code to see if it generates the same CRC from the
    % SIGNAL, SERVICE, LENGTHEXT bits acquired above
    % receivedCRC = headerRx(33:end);
    input = [signal, service, lengthext];
    % 1 + x^5 + x^12 + x^16
    poly = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 1]; % same polynomial as in Tx
    init = ones([1,15]);
    h = crc.generator('Polynomial',poly, ...
                        'InitialState',init);
    crcEncoded = generate(h, input');
    crcOnesComp = int8(not(num2str(crcEncoded) - '0')); % one's complement
    crcBits = fliplr(crcOnesComp(end-15:end)');

    if headerTx(33:end) ~= crcBits % headerRx(33:end)
        error('crc not matching!')
%     else
%         disp('crc checked')
    end
    
    % if it doesn't pass all error checking, return false
    if ~valid; checked = false; end
end
    