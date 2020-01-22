        % Calculate Median SMP Value for RammSonde Depth
%         tmpStrength = [];
        for ii = smpIx
            indx = find(invSMP{ii}.z <= depth);
            indx2 = find(SMP{ii}.depth_F <= depth);
            for jj = 1:11
                medM{nn}(ii,jj,1) = quantile(invSMP{ii}.M(indx,jj,2),q);
                medM{nn}(ii,jj,2) = quantile(invSMP{ii}.M(indx,jj,2),0.05);
                %                 medM{nn}(ii,jj,3) = quantile(invSMP{ii}.M(indx,jj,2),0.34);
                %                 medM{nn}(ii,jj,4) = quantile(invSMP{ii}.M(indx,jj,2),0.68);
                medM{nn}(ii,jj,3) = quantile(invSMP{ii}.M(indx,jj,2),0.95);
%                 if jj == 7 % Get Strength Data
%                     Sbox{nn} = [tmpStrength; invSMP{ii}.M(indx,jj,2)];
%                 end
                % Normalize SMP measures
%                 normMedM{nn}(ii,jj) = medM{nn}(ii,jj,1)./median(medM{nn}(ii,jj,1));
%                 subNormMedM{nn}(ii,jj) = medM{nn}(ii,jj,1) - median(medM{nn}(ii,jj,1));
            end
            for jj = 1:5
                medM{nn}(ii,jj+11,1) = quantile(invSMP{ii}.M2(indx,jj),q);
                medM{nn}(ii,jj+11,2) = quantile(invSMP{ii}.M2(indx,jj),0.05);
                %                 medM{nn}(ii,jj+11,3) = quantile(invSMP{ii}.M2(indx,jj),0.34);
                %                 medM{nn}(ii,jj+11,4) = quantile(invSMP{ii}.M2(indx,jj),0.68);
                medM{nn}(ii,jj+11,3) = quantile(invSMP{ii}.M2(indx,jj),0.95);
                % Normalize SMP measures
%                 normMedM{nn}(ii,jj+11) = medM{nn}(ii,jj+11,1)./median(medM{nn}(ii,jj+11,1));
%                 subNormMedM{nn}(ii,jj+11) = medM{nn}(ii,jj+11,1) - median(medM{nn}(ii,jj+11,1));
            end
            medM{nn}(ii,17,1) = quantile(SMP{ii}.force(indx2),q);
            medM{nn}(ii,17,2) = quantile(SMP{ii}.force(indx2),0.05);
            %             medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx),0.34);
            %             medM{nn}(ii,17,4) = quantile(SMP{ii}.force(indx),0.68);
            medM{nn}(ii,17,3) = quantile(SMP{ii}.force(indx2),0.95);
                % Normalize SMP measures
%                 normMedM{nn}(ii,17) = medM{nn}(ii,17,1)./median(medM{nn}(ii,17,1));
%                 subNormMedM{nn}(ii,17) = medM{nn}(ii,17,1) - median(medM{nn}(ii,17,1));            
            
        end
        zedIx = find(any(medM{nn}(:,:,1),2));
        dum = zeros(length(zedIx),17,3);
        for kk = 1:size(medM{nn},3)
            tmp = medM{nn}(zedIx,:,kk);
            dum(:,:,kk) = tmp;
        end
        medM{nn} = dum;
        clear dum
        if  ismember(nn,[3:3:24])
for jj = 1:17 
    %try 2
    tmp = [medM{nn-2}(:,jj,1);medM{nn-1}(:,jj,1);medM{nn}(:,jj,1)];
    siteMedian = nanmedian(tmp);
    for kk = nn-2:1:nn
    normMedM{kk}(:,jj) = medM{kk}(:,jj,1)./siteMedian;
    subNormMedM{kk}(:,jj) = medM{kk}(:,jj,1) - siteMedian;
    end
%try 1 unsuccessful
% normMedM{nn}(:,jj) = medM{nn}(:,jj,1)./median(medM{nn}(:,jj,1));
% subNormMedM{nn}(:,jj) = medM{nn}(:,jj,1) - median(medM{nn}(:,jj,1));
end
        end
        