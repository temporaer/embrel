function v = uptri_to_x(mat)

v=[];
for i=1:length(mat)
    v = [v mat(i,i:end)];
end
