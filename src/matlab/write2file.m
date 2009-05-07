% vim:ft=matlab
function write2file(PHIX, PSIY, COL)
  %MOL = PHIX{2};
  %CLA = PHIX{1};
  CLA = PHIX;
  FEA = [PSIY COL];
  %save /tmp/erl/mol.txt MOL -double -ascii
  save /home/schulzha/checkout/embrel/build/release/cla.txt CLA -double -ascii
  save /home/schulzha/checkout/embrel/build/release/fea.txt FEA -double -ascii
end

