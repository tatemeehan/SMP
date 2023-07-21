function [qc] = qcSMP(files)
% qcSMP is an interactive data review process
nfiles = length(files);
for kk = 1:nfiles
    figure(kk)
    msg = "Quality Control SMP Signal";
    opts = ["Good" "Drift" "Dry Run" "Bad"];
    file(kk,:) = string(files(kk).name);
%     choice(kk) = menu(msg,opts);
%     quality(kk,:) = opts(choice(kk));
    quality(kk,:) = opts(1);


end

qc = table(file,quality,'VariableNames',{'File Name', 'Quality'});