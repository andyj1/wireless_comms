classdef Filter
    properties(Constant)
        h = [];
    end
    
    methods(Static)
        
        function [h, upsampledChips, chipFilterDelay] = PulseShapeFilter(symbols, spc)
            span = 5;
            beta = 0.7;
            fcutoff = 7e6;
            fsamp = 88e6;
            M = span * spc; % filter order
            
            delayFilter = M/2; % default delay for this root raised cosine
            KaiserWindow = kaiser(M+1,1);
            % the most recent version of 'firrcos' is 'rcosdesign,' but
            % this method doesn't support truncating with a particular 
            % window other than 'rectangular', while we want to
            % use a Kaiser window as it is known to work rather well. So
            h = firrcos(M,fcutoff,beta,fsamp, ...
                        'rolloff','sqrt',delayFilter,KaiserWindow);
            
            % plot impulse response
            % fvtool(h, fsamp, 'Analysis', 'impulse');
            
            % compute filter delay 
            Ntap = length(h);
            sampleFilterDelay = Ntap-1;
            chipFilterDelay = sampleFilterDelay/spc;

            % upsample, normalize power by sqrt of sampling rate
            upsampledChips = zeros([spc*length(symbols),1]);
            for n = 0:1:size(symbols,1)-1
                upsampledChips(spc*n+1,1) = symbols(n+1,1)/spc;
            end
        end
        
        function [downsampledChips,bitDelay] = ...
            Receiver(filtTxSig, spc, chipFilterDelay, bps, chipSpreadLength)
            
            % take downsampled chips by spc rate
            chips = filtTxSig(1:spc:end, 1);
            
            chipShift = 0;
            % compute (if misaligned) alignment delay
            if rem(chipFilterDelay, chipSpreadLength) > 0
                chipShift = chipSpreadLength - rem(chipFilterDelay, chipSpreadLength);
            end
            
            % chipFilterDelay: from Tx filter
            % chipShift: from misalignment
            symbolDelay = (chipFilterDelay + chipShift) / chipSpreadLength;
            bitDelay = symbolDelay * bps;
            
            % align chips, make a column vector
            downsampledChips = [zeros([chipShift,1]);
                                chips(1:end-chipShift, 1)];
        end
    end % class methods
end % end classdef
