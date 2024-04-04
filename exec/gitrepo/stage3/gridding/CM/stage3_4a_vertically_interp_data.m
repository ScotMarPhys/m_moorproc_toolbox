function RT_merg_CM = stage3_4a_vertically_interp_data(Ufss,Vfss,Wfss,Pfss,pgg,JG)
    maxdepth = max(nanmean(Pfss,2));
    inodata = find(pgg>maxdepth);

    UGfs = nan(length(pgg),length(JG));
    VGfs = nan(length(pgg),length(JG));
    WGfs = nan(length(pgg),length(JG)); 
    UGfs_akima = nan(length(pgg),length(JG));
    VGfs_akima = nan(length(pgg),length(JG));
    WGfs_akima = nan(length(pgg),length(JG));   
    UGfs_pchip = nan(length(pgg),length(JG));
    VGfs_pchip = nan(length(pgg),length(JG));
    WGfs_pchip = nan(length(pgg),length(JG));          
    for ijj=1:length(JG)
        iuok = find(~isnan(Ufss(:,ijj)));
        ivok = find(~isnan(Vfss(:,ijj)));
        iwok = find(~isnan(Wfss(:,ijj)));
        if length(iuok)>1
            UGfs(:,ijj) = interp1(Pfss(iuok,ijj),Ufss(iuok,ijj),pgg,'linear') ;
            UGfs(inodata,ijj) = nan;
            UGfs_pchip(:,ijj) = interp1(Pfss(iuok,ijj),Ufss(iuok,ijj),pgg,'pchip') ;
            UGfs_pchip(inodata,ijj) = nan;
            if length(iuok)>2
                UGfs_akima(:,ijj) = akima(Pfss(iuok,ijj),Ufss(iuok,ijj),pgg) ;
                UGfs_akima(inodata,ijj) = nan;
            end
        end
        if length(ivok)>1
            VGfs(:,ijj) = interp1(Pfss(ivok,ijj),Vfss(ivok,ijj),pgg,'linear') ;    
            VGfs(inodata,ijj) = nan;        
            VGfs_pchip(:,ijj) = interp1(Pfss(ivok,ijj),Vfss(ivok,ijj),pgg,'pchip') ;  
            VGfs_pchip(inodata,ijj) = nan;
            if length(ivok)>2
                VGfs_akima(:,ijj) = akima(Pfss(iuok,ijj),Vfss(iuok,ijj),pgg) ;
                VGfs_akima(inodata,ijj) = nan;
            end
        end
        if length(iwok)>1
            WGfs(:,ijj) = interp1(Pfss(iwok,ijj),Wfss(iwok,ijj),pgg,'linear') ;    
            WGfs(inodata,ijj) = nan;
            WGfs_pchip(:,ijj) = interp1(Pfss(iwok,ijj),Wfss(iwok,ijj),pgg,'pchip') ;   
            WGfs_pchip(inodata,ijj) = nan;
            if length(iwok)>2
                WGfs_akima(:,ijj) = akima(Pfss(iuok,ijj),Wfss(iuok,ijj),pgg) ;
                WGfs_akima(inodata,ijj) = nan;
            end
        end  
    end
    % Alloacte variables into a structure
    RT_merg_CM.JG     = JG;
    RT_merg_CM.Ufs    = Ufs;
    RT_merg_CM.Vfs    = Vfs;
    RT_merg_CM.Wfs    = Wfs;
    RT_merg_CM.Pfs    = Pfs;    
    RT_merg_CM.PGfs   = pgg';    
    RT_merg_CM.UGfs   = UGfs;
    RT_merg_CM.VGfs   = VGfs;
    RT_merg_CM.WGfs   = WGfs;
    RT_merg_CM.UGfs_akima   = UGfs_akima;
    RT_merg_CM.VGfs_akima   = VGfs_akima;
    RT_merg_CM.WGfs_akima   = WGfs_akima;
    RT_merg_CM.UGfs_pchip   = UGfs_pchip;
    RT_merg_CM.VGfs_pchip   = VGfs_pchip;
    RT_merg_CM.WGfs_pchip   = WGfs_pchip;
    
    RT_merg_CM.comment{1,1}= 'JG -- julian day';
    RT_merg_CM.comment{2,1}= 'Ufs -- original stacked zonal velocity data from the deployments';
    RT_merg_CM.comment{3,1}= 'Vfs -- original stacked meridional data from the deployments ';
    RT_merg_CM.comment{4,1}= 'Wfs -- original stacked vertical data from the deployments ';    
    RT_merg_CM.comment{5,1}= 'Pfs -- original stacked Pressure data from the deployments';    
    RT_merg_CM.comment{6,1}= 'PGfs -- pressure grid ';
    RT_merg_CM.comment{7,1}= 'UGfs -- zonal velocity interpolated onto the pressure grid (PGfs) with a linear ';
    RT_merg_CM.comment{8,1}= 'VGfs -- meridional velocity interpolated onto the pressure grid (PGfs) with a linear ';
    RT_merg_CM.comment{9,1}= 'WGfs -- vertical velocity interpolated onto the pressure grid (PGfs) with a linear ';    
    RT_merg_CM.comment{10,1}= 'UGfs_akima -- zonal velocity interpolated onto the pressure grid (PGfs) with akima method';
    RT_merg_CM.comment{11,1}= 'VGfs_akima -- meridional velocity interpolated onto the pressure grid (PGfs) with akima method';
    RT_merg_CM.comment{12,1}= 'WGfs_akima -- vertical velocity interpolated onto the pressure grid (PGfs) with akima method';    
    RT_merg_CM.comment{13,1}= 'UGfs_pchip -- zonal velocity interpolated onto the pressure grid (PGfs) with a pchip ';
    RT_merg_CM.comment{14,1}= 'VGfs_pchip -- meridional velocity interpolated onto the pressure grid (PGfs) with a pchip ';
    RT_merg_CM.comment{15,1}= 'WGfs_pchip -- vertical velocity interpolated onto the pressure grid (PGfs) with a pchip '; 
end