from netCDF4 import Dataset

def output_init(self):
    fh01 = Dataset(self.output, mode='w',format="NETCDF4")
    fh01.createDimension('z_r',self.nz)
    fh01.createDimension('z_w',self.nz+1)
    fh01.createDimension('time',None)
    fh01.fcorSIN         = self.fcor
    fh01.fcorCOS         = self.ecor
    fh01.eddy_diff       = str(self.eddy_diff)
    fh01.evd             = str(self.ED_evd)
    fh01.mass_flux_tra   = str(self.MF_tra)
    fh01.mass_flux_dyn   = str(self.MF_dyn)
    fh01.mass_flux_tke   = str(self.MF_Etransfer_to_tke)
    fh01.mass_flux_tke_trplCorr = str(self.MF_tke_trplCorr)
    fh01.cp              = self.cp
    fh01.rho0            = self.eos_params[0]
    fh01.alpha           = self.eos_params[1]
    fh01.beta            = self.eos_params[2]
    fh01.T0              = self.eos_params[3]
    fh01.S0              = self.eos_params[4]
    #self.min_Threshold
    fh01.tkemin          = self.min_Threshold[0]
    fh01.akvmin          = self.min_Threshold[1]
    fh01.aktmin          = self.min_Threshold[2]
    fh01.mxlmin          = self.min_Threshold[3]
    # 
    ocean_time = fh01.createVariable('ocean_time','f8',('time')); ocean_time[:] = 0.
    ocean_time.units = 'seconds'
    ocean_time.long_name = 'time since initialization'
    taux = fh01.createVariable('taux','f8',('time')); taux[:] = self.ustr_sfc
    tauy = fh01.createVariable('tauy','f8',('time')); tauy[:] = self.vstr_sfc
    Qns  = fh01.createVariable('Qns','f8',('time'));   Qns[:] = self.stflx[self.itemp]*self.rho0*self.cp           #non-solar
    Qs   = fh01.createVariable( 'Qs','f8',('time'));    Qs[:] = self.srflx*self.rho0*self.cp                          #solar
    Fw   = fh01.createVariable( 'Fw','f8',('time'));    Fw[:] = self.stflx[self.isalt]            #freshwater
    zw   = fh01.createVariable('z_w','f8',('z_w')); zw[:] = self.z_w[:]
    zr   = fh01.createVariable('z_r','f8',('z_r')); zr[:] = self.z_r[:]
    var  = fh01.createVariable('u','f8',('time','z_r')); var[0,:] = self.u_n[:]; var.units = 'm s-1'; del var
    var  = fh01.createVariable('v','f8',('time','z_r')); var[0,:] = self.v_n[:]; var.units = 'm s-1'; del var
    var  = fh01.createVariable('temp','f8',('time','z_r')); var[0,:] = self.t_n[:,self.itemp]; var.units = 'Celsius'; del var
    var  = fh01.createVariable('salt','f8',('time','z_r')); var[0,:] = self.t_n[:,self.isalt]; var.units = 'psu'; del var
    var  = fh01.createVariable('Etot','f8',('time')); var[0] = self.vint_Etot
    var  = fh01.createVariable('Ekin','f8',('time')); var[0] = self.vint_Ekin
    var  = fh01.createVariable('Epot','f8',('time')); var[0] = self.vint_Epot
    var  = fh01.createVariable('Etke','f8',('time')); var[0] = self.vint_TKE
    var  = fh01.createVariable('Eeps','f8',('time')); var[0] = self.vint_Eps
    var  = fh01.createVariable('hmxl','f8',('time')); var[0] = self.hmxl
    var  = fh01.createVariable('WT','f8',('time','z_w')); var[0,:] = self.wted[:]+self.wtmf[:]; var.units = 'K m s-1'; del var
    var  = fh01.createVariable('WU','f8',('time','z_w')); var[0,:] = self.wued[:]+self.wumf[:]; var.units = 'm2 s-2'; del var
    var  = fh01.createVariable('WV','f8',('time','z_w')); var[0,:] = self.wved[:]+self.wvmf[:]; var.units = 'm2 s-2'; del var
    var  = fh01.createVariable('WTKE','f8',('time','z_r')); var[0,:] = self.wtke[:] ; var.units = 'm3 s-3'; del var
    #
    if self.eddy_diff:
        var  = fh01.createVariable('tke','f8',('time','z_w')); var[0,:] = self.tke_n[:]; var.units = 'm2 s-2'; del var
        var  = fh01.createVariable('Akv','f8',('time','z_w')); var[0,:] = self.akv[:]; var.units = 'm2 s-1'; del var
        var  = fh01.createVariable('Akt','f8',('time','z_w')); var[0,:] = self.akt[:]; var.units = 'm2 s-1'; del var
        var  = fh01.createVariable('bvf','f8',('time','z_w')); var[0,:] = self.bvf[:]; var.units = 's-2'; del var
        var  = fh01.createVariable('lup','f8',('time','z_w')); var[0,:] = self.lupw[:]; var.units = 'm'; del var
        var  = fh01.createVariable('ldw','f8',('time','z_w')); var[0,:] = self.ldwn[:]; var.units = 'm'; del var
        var  = fh01.createVariable('WT_ED','f8',('time','z_w')); var[0,:] = self.wted[:]; var.units = 'K m s-1'; del var
        var  = fh01.createVariable('WU_ED','f8',('time','z_w')); var[0,:] = self.wued[:]; var.units = 'm2 s-2'; del var
        var  = fh01.createVariable('WV_ED','f8',('time','z_w')); var[0,:] = self.wved[:]; var.units = 'm2 s-2'; del var
        var  = fh01.createVariable('Prdtl','f8',('time','z_w')); var[0,:] = self.akv[:]/self.akt[:]; var.units = 'none'; del var
        var  = fh01.createVariable('epsil','f8',('time','z_w')); var[0,:] = self.eps_n[:]; var.units = 'm2 s-3'; del var
        var  = fh01.createVariable('Buoy_prod','f8',('time','z_w')); var[0,:] = self.Bprod[:]; var.units = 'm2 s-3'; del var
    if self.MF_tra:
        var = fh01.createVariable('a_p','f8',('time','z_w')); var[0,:] = self.ap[:]; del var
        var = fh01.createVariable('zinv','f8',('time')); var[0] = self.zinv
        var = fh01.createVariable('w_p','f8',('time','z_w')); var[0,:] = self.wp[:]; del var
        var = fh01.createVariable('B_p','f8',('time','z_w')); var[0,:] = self.buoyMF[:]; del var
        var = fh01.createVariable('temp_p','f8',('time','z_w')); var[0,:] = self.tp[:,self.itemp]; del var
        var = fh01.createVariable('salt_p','f8',('time','z_w')); var[0,:] = self.tp[:,self.isalt]; del var
        var = fh01.createVariable('Ent','f8',('time','z_r')); var[0,:] = self.ent[:]; del var
        var = fh01.createVariable('Det','f8',('time','z_r')); var[0,:] = self.det[:]; del var
        var  = fh01.createVariable('WT_MF','f8',('time','z_w')); var[0,:] = self.wtmf[:]; del var
        var  = fh01.createVariable('WU_MF','f8',('time','z_w')); var[0,:] = self.wumf[:]; del var
        var  = fh01.createVariable('WV_MF','f8',('time','z_w')); var[0,:] = self.wvmf[:]; del var
        var = fh01.createVariable('u_p','f8',('time','z_w')); var[0,:] = self.up[:]; del var
        var = fh01.createVariable('v_p','f8',('time','z_w')); var[0,:] = self.vp[:]; del var
        if self.MF_Etransfer_to_tke:
            var = fh01.createVariable('buoyMF','f8',('time','z_w')); var[0,:] = self.buoyMF[:]; del var
            var = fh01.createVariable('shearMF','f8',('time','z_w')); var[0,:] = self.shearMF[:]; del var
        if self.MF_tke_trplCorr:
#                var = fh01.createVariable('d_we_dz_MF','f8',('time','z_w')); var[0,:] = self.triple_corr[:]; del var
            var = fh01.createVariable('tke_p','f8',('time','z_w')); var[0,:] = self.tkep[:]; del var
        var = fh01.createVariable('vort_p','f8',('time','z_w')); var[0,:] = self.vortp[:]; del var
    fh01.close()
#












#
def output_state(self,TimeInSecs,kout):
    fh01 = Dataset(self.output, mode='a',format="NETCDF4")
    fh01.variables['ocean_time'][kout] = TimeInSecs
    fh01.variables['taux'][kout]       = self.ustr_sfc
    fh01.variables['tauy'][kout]       = self.vstr_sfc
    fh01.variables['Qns'][kout]        = self.stflx[self.itemp]*self.cp*self.rho0
    fh01.variables['Qs'][kout]         = self.srflx
    fh01.variables['Fw'][kout]         = self.stflx[self.isalt]*self.cp*self.rho0
    fh01.variables['u'][kout,:]        = self.u_n[:]
    fh01.variables['v'][kout,:]        = self.v_n[:]
    fh01.variables['temp'][kout,:]     = self.t_n[:,self.itemp]
    fh01.variables['salt'][kout,:]     = self.t_n[:,self.isalt]
    fh01.variables['Etot'][kout]       = self.vint_Etot
    fh01.variables['Ekin'][kout]       = self.vint_Ekin
    fh01.variables['Epot'][kout]       = self.vint_Epot
    fh01.variables['Etke'][kout]       = self.vint_TKE
    fh01.variables['Eeps'][kout]       = self.vint_Eps
    fh01.variables['hmxl'][kout]       = self.hmxl
    fh01.variables['WT'][kout,:]       = self.wted[:]+self.wtmf[:]
    fh01.variables['WU'][kout,:]       = self.wued[:]+self.wumf[:]
    fh01.variables['WV'][kout,:]       = self.wved[:]+self.wvmf[:]
    fh01.variables['Buoy_prod'][kout,:]= self.Bprod[:]
    fh01.variables['WTKE'][kout,:]= self.wtke[:]

    if self.eddy_diff:
        fh01.variables['tke'][kout,:] = self.tke_n[:]
        fh01.variables['Akv'][kout,:] = self.akv[:]
        fh01.variables['Akt'][kout,:] = self.akt[:]
        fh01.variables['bvf'][kout,:] = self.bvf[:]
        fh01.variables['lup'][kout,:] = self.lupw[:]
        fh01.variables['ldw'][kout,:] = self.ldwn[:]
        fh01.variables['WT_ED'][kout,:] = self.wted[:]
        fh01.variables['WU_ED'][kout,:] = self.wued[:]
        fh01.variables['WV_ED'][kout,:] = self.wved[:]
        fh01.variables['Prdtl'][kout,:] = self.akv[:]/self.akt[:]
        fh01.variables['epsil'][kout,:] = self.eps_n[:]
    if self.MF_tra:
        fh01.variables['zinv'][kout]  = self.zinv
        fh01.variables['a_p'][kout,:] = self.ap[:]
        fh01.variables['w_p'][kout,:] = self.wp[:]
        fh01.variables['B_p'][kout,:] = self.Bp[:]
        fh01.variables['temp_p'][kout,:] = self.tp[:,self.itemp]
        fh01.variables['salt_p'][kout,:] = self.tp[:,self.isalt]
        fh01.variables['Ent'][kout,:] = self.ent[:]
        fh01.variables['Det'][kout,:] = self.det[:]
        fh01.variables['WT_MF'][kout,:] = self.wtmf[:]
        fh01.variables['WU_MF'][kout,:] = self.wumf[:]
        fh01.variables['WV_MF'][kout,:] = self.wvmf[:]
        fh01.variables['u_p'][kout,:] = self.up[:]
        fh01.variables['v_p'][kout,:] = self.vp[:]
        if self.MF_Etransfer_to_tke:
            fh01.variables['buoyMF'][kout,:] = self.buoyMF[:]
            fh01.variables['shearMF'][kout,:] = self.shearMF[:]
            fh01.variables['Buoy_prod'][kout,:]= self.Bprod[:] + self.buoyMF[:] #ED+MF buoyancy production
        if self.MF_tke_trplCorr:
#                fh01.variables['d_we_dz_MF'][kout,:] = self.triple_corr[:]
            fh01.variables['tke_p'][kout,:] = self.tkep[:]
        fh01.variables['vort_p'][kout,:] = self.vortp[:]
    fh01.close()

