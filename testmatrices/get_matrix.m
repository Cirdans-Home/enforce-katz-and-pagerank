clear all
close all 
clc

pb_class    = {'HB/'};
pb_name   = {'gre_1107'};

for k = 1:length(pb_name)  
Problem     = ssget([pb_class{k},pb_name{k}]);
save(['./',pb_name{k}], "Problem");
end