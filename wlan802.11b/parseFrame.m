function [preamble, header] = parseFrame(RxBits)    
    pos = 1;
    preamble = RxBits(pos:144);
    pos = pos + 144;
    header = RxBits(pos:pos+48-1);
end