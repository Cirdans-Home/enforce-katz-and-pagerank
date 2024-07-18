clear all
close all 
clc

pb_class    = {'vanHeukelum/','AG-Monien/'};
pb_name   = {'cage10','se'};

for k = 1:length(pb_name)  
Problem     = ssget([pb_class{k},pb_name{k}]);
save(['./',pb_name{k}], "Problem");
end