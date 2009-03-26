% vim:ft=matlab
1;
cd /tmp;
%feat_feat_mat;
%feat_mol_mat;
%feat_klass_mat;
code_data;

feat_klass = feat_klass / sum(feat_klass(:)); % make probabilities
feat_feat  = feat_feat / sum(feat_feat(:));   % make probabilities
%feat_mol   = feat_mol / sum(sum(feat_mol));   % make probabilities


cfg = '/home/sarx/prog/uni/embrel/src/matlab/config.cfg';

cd /home/sarx/prog/uni/embrel/src/third_party/coemb
%o = struct('n_restarts', 10, 'x_cond', 1, 'x_marg', 'M', 'y_marg', 'M', 'b_fix_psi', 0);
%o = struct('n_restarts', 2, 'b_cond', 0, 'x_marg', 'U', 'y_marg', 'U', 'b_fix_psi', 0, 'w_nxy', [1,0.30]);
o = struct();
o.n_restarts = str2num(ml_GetPrivateProfileString('CODE','n_restarts', cfg));
o.a_cond     = 0;
o.b_cond     = str2num(ml_GetPrivateProfileString('CODE','b_cond'    , cfg)); % 1
o.x_marg     = ml_GetPrivateProfileString('CODE','x_marg'    , cfg); % M
o.y_marg     = ml_GetPrivateProfileString('CODE','y_marg'    , cfg); % F;
if(str2num(ml_GetPrivateProfileString('CODE','use_pxx', cfg)) == 1)
	o.pxx0       = feat_feat;
	%o.w_pxx0 = 0.5*size(feat_klass,2)/size(feat_klass,1);
	%o.w_pxx0 = 10.0;
	o.w_pxx0 = eval(ml_GetPrivateProfileString('CODE', 'w_pxx0', cfg));
end

%[PHIX, PSIY, L] = code({pxy_data'}, 2, o);
%phi_x = PHIX{2};
%psi_y = PSIY{2};
%[PHIX, PSIY, L] = code({feat_klass, feat_mol}, 2, o);
[PHI, PSI, L] = code({feat_klass}, 2, o);
psi = PHI{2};
phi = PSI{2};

clf;
hold on;

s = sum(feat_klass);
%pxycolor = max(pxy_data./[s;s;s;s;s],[],1);
pxycolor = max(feat_klass',[],1);
hold off;
cd /home/sarx/prog/uni/ma/molemb1/matlab
write2file(phi, psi, pxycolor');
