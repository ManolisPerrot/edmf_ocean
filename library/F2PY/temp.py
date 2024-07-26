        if self.bc_P09 == 'consistent': MF_sfc_flux = 0. #bizarre?

            #===================================================
            # Advance tracers to n+1 (vertical diffusion only)
            #===================================================
            #
            stflx0[self.itemp] = self.stflx[self.itemp]
            if self.bc_P09 =='consistent':
                #print("=================")
                #print(self.stflx[self.itemp])
                stflx0[self.itemp] = stflx0[self.itemp] - MF_sfc_flux

                                if self.bc_P09 == 'consistent': MF_sfc_flux = self.Fmass[-1]*(self.tp[-1,self.itemp]-self.t_np1_star[-1,self.itemp])
                
comprend pas à quoi sert compute_MF_forcing                