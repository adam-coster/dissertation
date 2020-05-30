% Measuring unbiased mutual information beween ligands and readouts for
% TGFB3 and BMP4

ligands = {'TGFB3','TGFB1','BMP4','Wnt3A'};
features = {'median','sum'} ;

savename = 'ligandInformation.csv' ;
f = fopen(savename, 'w') ;
for rIdx = 1:length(ligands),
    for fIdx = 1:length(features),
        fprintf(f, [ligands{rIdx} '.' features{fIdx} '_nucleus.meanMI,' ...
                    ligands{rIdx} '.' features{fIdx} '_nucleus.sdMI,'] ) ;
    end
end
fprintf(f,'\n') ;

for rIdx = 1:length(ligands),
    for fIdx = 1:length(features),
        files = dir(['.' ligands{rIdx} '*.' features{fIdx} '_nucleus.csv']) ;
        % Grab the file contents (skipping the header line)
        mi = mutual_information( {csvread(files(1).name,2)',csvread(files(2).name,2)'},...
                                 10:5:60,[5:.5:10]/10 ) ;
        fprintf(f,[num2str(mean(mi)) ',' num2str(std(mi)) ',' ]) ;
    end
end
fprintf(f,'\n')
fclose(f)
