function [checked] = checkFrame(preamble, header)
    
%     header'
    % check preamble
    scramInit = [1 1 0 1 1 0 0]';
    syncDescrambled = wlanScramble(preamble(1:128),scramInit);
    if syncDescrambled ~= ones([128,1])
%         error('sync not matching!')
    end
    
    sfd = preamble(129:129+16-1);
    if sfd ~= fliplr([1 1 1 1 0 0 1 1 1 0 1 0 0 0 0 0])
%         error('sfd not matching!')
    end
    
    % check header
    signal = header(1:8);
    if ~isequal(signal, [zeros(1,4), 1 0 1 0]) || ...
        ~isequal(signal, [zeros(1,3), 1 0 1 0 0]) || ...
        ~isequal(signal, [zeros(1,2), 1 1 0 1 1 1]) || ...
        ~isequal(signal, [zeros(1,1), 1 1 0 1 1 1 0])
      
%         error('signal not matching!')
    end
    
    service = header(9:9+8-1);
    % check modulation scheme (CCK)
    if service(4) ~= 0
%         error('MCS is not CCK!')
    end
    
    length = header(17:17+16-1);
    
    crccode = header(33:end);
    input = [signal, service, length];
    % 1 + x^5 + x^12 + x^16
    h = crc.generator('Polynomial',[1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 1], ...
                        'InitialState',ones([1,15]));
    crcEncoded = generate(h, input');
    crcOnesComp = int8(not(num2str(crcEncoded) - '0'));
    crcBits = fliplr(crcOnesComp(end-15:end)');
    if crccode ~= crcBits
%         error('crc not matching!')
    end
    checked = true;    
end
    