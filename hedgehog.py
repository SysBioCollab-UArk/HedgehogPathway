from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

Model()

# Monomers
# Parameters
# Initials
# Rules
# Observables
# simulation commands

# Absence of Hh Signaling Pathway
#PTCH is bind to the primary cilia (for now)
#Name is capitalized, and lower case whatever binds
Monomer('PTCH',['hh'])
Monomer('HH', ['ptch'])
Monomer('SUFU', ['gli'])
Monomer('Gli',['sufu','state'],{'state' : ['U', 'P','PP','PPP']})
#Do I also need to represent the phosphorylation of Gli due to its interaction with other enzymes?
Monomer('CK1')
Monomer('PKA')
Monomer('GSK3B')
Monomer('GliR',['b'])
Monomer('GliA')


Parameter('kf_ptch_hh',1)
Parameter('kr_ptch_hh',1)
Parameter('kf_gli_sufu', 1)
Parameter('kr_gli_sufu', 1)
Parameter('k_gli_phos_ck1', 1)
Parameter('k_gli_phos_pka', 1)
Parameter('k_gli_phos_gsk3b', 1)
Parameter('k_gli_gliR', 1)

#where do I represent the conversion of Gli into GliA and GliR?(Parameter section?)
Rule('PTCH_binds_HH',PTCH(hh=None) + HH(ptch=None) | PTCH(hh=1) % HH(ptch=1), kf_ptch_hh, kr_ptch_hh)
Rule( 'Gli_binds_SUFU', Gli(sufu=None,state='U') + SUFU(gli=None) | Gli(sufu=1, state='U') % SUFU(gli=1), kf_gli_sufu, kr_gli_sufu)
Rule('Gli_phos_CK1', Gli(sufu=1, state='U') % SUFU(gli=1) + CK1() >> Gli(sufu=1, state='P') % SUFU(gli=1) + CK1() , k_gli_phos_ck1)
Rule('Gli_phos_PKA', Gli(sufu=1, state='P') % SUFU(gli=1) + PKA() >> Gli(sufu=1, state='PP') % SUFU(gli=1) + PKA() , k_gli_phos_pka)
Rule('Gli_phos_GSK3B', Gli(sufu=1, state='PP') % SUFU(gli=1) + GSK3B() >> Gli(sufu=1, state='PPP') % SUFU(gli=1) + GSK3B() , k_gli_phos_gsk3b)
Rule('Gli_To_GliR', Gli(sufu=1, state='PPP') % SUFU(gli=1) >> GliR() + SUFU(gli=None), k_gli_gliR)

#Parameters for initials
Parameter('PTCH_init', 1)
Parameter('HH_init', 1)
Parameter('Gli_init', 1)
Parameter('SUFU_init', 1)
Parameter('CK1_init', 1)
Parameter('PKA_init', 1)
Parameter('GSK3B_init', 1)

#Initials
Initial(PTCH(hh=None), PTCH_init)
Initial(HH(ptch=None), HH_init)
Initial(Gli(sufu=None, state='U'), Gli_init)
Initial(SUFU(gli=None), SUFU_init)
Initial(CK1(), CK1_init)
Initial(PKA(), PKA_init)
Initial(GSK3B(), GSK3B_init)


# Observable('A_free', A(b=None))
# Observable('A_phos', A(state='P'))
# # Observable('A_bound, A(b=1') % B(a=1))
# Observable('A_bound', A(b=ANY))

# Presence of Hh Signaling Pathway

