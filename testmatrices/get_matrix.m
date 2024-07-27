clear all
close all 
clc

pb_class    = {'JGD_SPG/','Arenas/','vanHeukelum/','vanHeukelum/','DIMACS10/',...
                      'DIMACS10/', 'DIMACS10/','DIMACS10/','DIMACS10/','DIMACS10/',...
                      'DIMACS10/','HB/','DIMACS10/','DIMACS10/','DIMACS10/'};
pb_name     = {'EX5','PGPgiantcompo','cage10','cage11','cs4','ct2010','cti','data',...
                        'de2010','delaunay_n16','fe_4elt2','gre_1107','nh2010','uk','vt2010'};


for k = 1:length(pb_name)  
Problem     = ssget([pb_class{k},pb_name{k}]);
save(['./',pb_name{k}], "Problem");
end