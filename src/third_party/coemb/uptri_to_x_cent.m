% Takes a matrix, and converts its upper triangle to 
% the set of variables we work with. It ignores the upper right corner of the matrix
% since this is not a variable, and is set so that all 
% the elements sum up to zero

function v = uptri_to_x_cent(mat)

v=[];
for i=1:length(mat)
    if i==1
        v = [v mat(i,i:end-1)];
    else
        v = [v mat(i,i:end)];
    end
end
