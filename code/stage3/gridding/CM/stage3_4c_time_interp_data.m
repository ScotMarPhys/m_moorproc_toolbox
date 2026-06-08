function RT_merg_CM = stage3_4c_time_interp_data(U,V,pgg,JG,idepth,RT_merg_CM)
    %--------
    % Apr 2024 - removed W as not sensible in merged product (K. Burmeister)
    %
    %--------
    [m,n] = size(RT_merg_CM.UGfs);
    UG_2 = NaN * ones(m,n); VG_2 = NaN * ones(m,n); 

    % ingore upper bins without data
    I = find(pgg == idepth);

    % copy the top idepth into the new file
    % no temporal interpolation
    UG_2([1:I-1],:)=U([1:I-1],:);
    VG_2([1:I-1],:)=V([1:I-1],:);


        for i = I: length(pgg) % for each depth

            % locate all non nan values in the despiked zonal velocity
            iu = find(~isnan(U(i,:)));
            if length(iu) < 2
                continue
            end
            % locate all non nan values in the despiked meridional velocity
            iv = find(~isnan(V(i,:)));

            % interpolate in time over the missing data
            VG_2(i,:) = interp1(JG(iv), V(i, iv), JG);
            % interpolate in time over the missing data
            UG_2(i,:) = interp1(JG(iu), U(i, iu), JG);

            % set NAN for nans gap of more than 10 days             
            ibad = find([diff(iv)]>20); 
            if length(ibad)>=1
                for isss=1:length(ibad)  
                    VG_2(i,iv(ibad(isss)):iv(ibad(isss)+1))=nan;
                    UG_2(i,iv(ibad(isss)):iv(ibad(isss)+1))=nan;
                end
            end
        end

    % KB, 27/3/2024
    % Added code to remove horizontal interpolation when mooring was knocked down 
    % This was not needed before because data was vertical extrapolated to surface 
    % first which is now done in later step.

    for i = 1:length(JG);
        iu = find(~isnan(RT_merg_CM.UGfs(:,i)),1,'first');
        iv = find(~isnan(RT_merg_CM.VGfs(:,i)),1,'first');

        VG_2(1:iv-1,i)=nan;
        UG_2(1:iu-1,i)=nan;
    end
    RT_merg_CM.UGfs2= UG_2;
    RT_merg_CM.VGfs2= VG_2;


    RT_merg_CM.comment{16,1}= ['UGfs2 -- zonal velocity interpolated ',...
                        'onto the time grid (JG) after despiking with a linear '];  
    RT_merg_CM.comment{17,1}= ['VGfs2 -- meridional velocity interpolated',...
                        ' onto the time grid (JG) after despiking with a linear '];

end