
n_bst=1000;
max_repeats=10;
dim_list = [2 3 4 7 10];
topics = [4 3];

i=1;str{i}='le'  ;func{i}='logemb_grad';prms{i}=max_repeats;
i=2;str{i}='ca'  ;func{i}='corr_anl'   ;prms{i}='';
i=3;str{i}='mdsc';func{i}='my_mds'     ;prms{i}='cosine';
i=4;str{i}='mdse';func{i}='my_mds'     ;prms{i}='euclide';
 

do = [0 1 0 0];
procmodes=[0:3];

for(i=1:4)  
  if(do(i))
    for(procmode=procmodes)      
      main_all(topics,dim_list,n_bst,str{i},func{i},prms{i},procmode);          
    end
  end
end
  