
i=1;LEN=1; clr{i}=[0.0 0.0 0.7];str{i}='le'; do{i}=1;mtr{i}='c';
i=2;CA=2; clr{i}=[0.7 0.0 0.0];str{i}='ca'; do{i}=1;mtr{i}='';
i=3;MDS=3;clr{i}=[0.0 0.7 0.5];str{i}='mds';do{i}=1;mtr{i}='ctf'
i=4;SVD=4;clr{i}=[0.5 0.0 0.5];str{i}='svd';do{i}=1;mtr{i}='';
i=5;MDS=5;clr{i}=[0.7 0.7 0.5];str{i}='mds';do{i}=1;mtr{i}='etf'
i=6;LEC=6; clr{i}=[0.0 0.0 0.7];str{i}='le'; do{i}=1;mtr{i}='c';


do_figures = [0 1 1 1 ];


LW = 'LineWidth';
FS = 'FontSize';
MS = 'MarkerSize';
CLR = 'Color';

% Set figure properties 
figure(1); clf; hold on ; 
set(gcf,'DefaultAxesFontSize',12);
set(gcf,'DefaultAxesFontName','helvetica');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0.5 0.6 26 6]);
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[1 19 0 0]);
y_low=0.2;y_len=0.7;
x_first=0.05;x_len=0.18;x_jump=0.25;

% Figure 1 - embedding
% ========
% Read pre calculated data
if(do_figures(1))
  nx=1000;
  infix='-1-2-12';
  filename = sprintf('CALE_topics%s_dims_nx_%d',infix,nx);
  load(filename,'le_psi');
  
  % Plot figure 
  subplot('position',[x_first y_low,x_len,y_len]); cla;hold on ;
  plot(le_psi(2001:3000,1),le_psi(2001:3000,2),'b.',MS,6);
  plot(le_psi(   1:1000,1),le_psi(   1:1000,2),'rx',MS,3);
  plot(le_psi(1001:2000,1),le_psi(1001:2000,2),'mo',MS,3);
  axis([0 1 0 1 ])
  box on ;
end

% Figure 2 - purity
% =================
% Read pre calculated data
if(do_figures(2))
  nx=1000;
  infix='-1-2-12';
  infix='-18-19';  
  subplot('position',...
	  [x_first+x_jump y_low,x_len,y_len]); cla; hold on ;
  for(i=1:4)
    if(do{i})
      filename = sprintf('CALE_%s%s_topics%s_dims_nx_%d',str{i},mtr{i},infix,nx)
      tmp=load(filename,sprintf('%s_purity',str{i}));          
      eval(sprintf('purity=tmp.%s_purity',str{i}));
      plot(purity{1}(1:1000),CLR,clr{i},LW,2);  
    end
  end
  
  % Plot figure 
  axis([0 1001 0.39 1])
  set(gca,'Xscale','log')
  set(gca,'Xtick',[1 10 100 1000])
  set(gca,'Xticklabel',[1 10 100 1000])
  set(gca,'Ytick',[0.4:0.1:1]);
  box on 
  hl=legend('CODE','CA','MDS','SVD'),set(hl,FS,8);
  xlabel('N nearest neighbors',FS,14);
  ylabel('purity',FS,14);
end


% Figure 3 - as a function of data size
% =====================================
if(do_figures(3))
  subplot('position',...
	  [x_first+x_jump*2 y_low,x_len,y_len]); cla;hold on ;
  
  i_dim=1;
  for(i=1:4)        
    nxs=[1000 2000 3000 4000];
    for(i_nx=1:length(nxs))
      nx= nxs(i_nx);      
      if(do{i})
	filename=sprintf('CALE_%s%s_topics%s_dims_nx_%d',...
			 str{i},mtr{i},infix,nx);
	tmp=load(filename,sprintf('%s_purity',str{i}));          
	eval(sprintf('purity=tmp.%s_purity',str{i}));	
	m(i_nx) = mean(purity{i_dim}(1:1000));
      end
    end
    plot(nxs,m,CLR,clr{i},LW,3);    
  
    set(gca,'Xtick',nxs)
    set(gca,'Ytick',[0.33 0.4:0.1:1]);
    box on 
    hl=legend('CODE','CA','MDS','SVD');set(hl,FS,8);
    xlabel('N words',FS,14);
    ylabel('purity',FS,14);
    axis([min(nxs)-200 max(nxs)+200 0.33 1]);
  end  
end
  
% Figure 4 - as a function of dimention
% =====================================
if(do_figures(4))
  nx=1000;
  dim_list = [2 3 4 5 6 8 10];
  infix='-1-2-12';
  
  subplot('position',...
	  [x_first+x_jump*3 y_low,x_len,y_len]); cla; hold on;
  
  for(i=1:4)      
    if(do{i})
      filename = sprintf('CALE_%s%s_topics%s_dims_nx_%d',...
			 str{i},mtr{i},infix,nx);
      tmp=load(filename,sprintf('%s_purity',str{i}));          
      eval(sprintf('purity=tmp.%s_purity',str{i}));
      for(i_dim=1:length(dim_list))
	m(i_dim)=mean(purity{i_dim}(1:1000));
      end
      plot(dim_list,m(1:length(dim_list)) ,CLR,clr{i}, LW,2);      
    end    
  end
  set(gca,'Xtick',dim_list)
  set(gca,'Ytick',[0.5:0.1:1]);
  hl=legend('CODE','CA','MDS','SVD');set(hl,FS,8);
  
  box on 
  xlabel('dimension',FS,14);
  ylabel('purity',FS,14);
  axis([1.8 10 0.6 1])
  
end


return

% Table
% =====
sets ={'-1-2-12' '-2-3' '-4-3' '-12-13' ,'-12-13-14','-12-13-14-15','-18-19'};

fprintf('%15s %4s %4s %4s\n', 'sets', 'CA', 'SVD','CODE');
for(i_set=1:length(sets))
  infix = sets{i_set};
  nx=1000; 
  filename = sprintf('CALE_ca_topics%s_dims_nx_%d',infix,nx);
  load(filename,'ca_purity');
  filename = sprintf('CALE_le_topics%s_dims_nx_%d',infix,nx);
  load(filename,'le_purity');  
  filename = sprintf('CALE_svd_topics%s_dims_nx_%d',infix,nx);
  load(filename,'svd_purity');  
  
  fprintf('%15s %4.2f %4.2f %4.2f\n', infix, ...
	  mean(ca_purity{1}(1:1000)), ...
	  mean(svd_purity{1}(1:1000)), ...	  
	  mean(le_purity{1}(1:1000)));
end

