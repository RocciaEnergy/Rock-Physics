# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 22:45:28 2022

@author: tamir
"""

import numpy as np

def CritPor(vp1,vs1,ro1,vp2,vs2,ro2,phicr):
    """
    Function computes velocities, elastic moduli, and bulk density at the 
    critical porosity.
    # =========================================================================
    #                           INPUTS
    # =========================================================================
    Takes the following inputs:
        vp1,vs1,vp2,vs2:   Velocities of the two constituents
        ro1,ro2:           Densities of the two constituents
        phicr:             Critical porosity
    # =========================================================================
    #                           OUTPUTS
    # =========================================================================
    Function returns a list of outputs:
        [vpcr,vscr,rocr,mcr,kcr,mucr]
        vpcr, vscr:       Velocities on critical porosity
        rocr:             Bulk density on critical porosity
        kcr, mucr:        Bulk and shear moduli on critical porosity
        rocr:             Bulk density on critical porosity
        mcr:              
    """

    m1=ro1*vp1**2; m2=ro2*vp2**2; mu1=ro1*vs1**2; mu2=ro2*vs2**2;
    k1=m1-(4/3)*mu1; k2=m2-(4/3)*mu2;
    
    mcr=(m1*m2)/((1-phicr)*m2+phicr*m1);
    mucr=(mu1*mu2)/((1-phicr)*mu2+phicr*mu1);
    kcr=(k1*k2)/((1-phicr)*k2+phicr*k1);
    rocr=(1-phicr)*ro1+phicr*ro2; vscr=np.sqrt(mucr/rocr);
    vpcr=np.sqrt((kcr+(4/3)*mucr)/rocr);
    
    outlistcr = [vpcr,vscr,rocr,mcr,kcr,mucr]
    
    #plt.plot(phicr,kcr,'b',phicr,mucr,'r')
    return outlistcr
