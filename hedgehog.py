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
Monomer('GliR',['b','loc'], {'loc' : ['cyt','nuc']})
Monomer('GliA',['b', 'state', 'loc'],{'state' : ['U', 'P'],'loc':[ 'cyt','nuc']})
Monomer('SMO',['pc'])



Parameter('kf_ptch_hh',1)
Parameter('kr_ptch_hh',1)
Parameter('kf_gli_sufu', 1)
Parameter('kr_gli_sufu', 1)
Parameter('k_gli_phos_ck1', 1)
Parameter('k_gli_phos_pka', 1)
Parameter('k_gli_phos_gsk3b', 1)
Parameter('k_gli_ppp_sufu', 1)
Parameter('k_gli_gliR', 1)
Parameter('k_gliA_deg', 1)
Parameter('k_gli_gliA_U',1)
Parameter('k_gli_gliA_P',1)
Parameter('k_gliRcyt_gliRnuc',1)
Parameter('k_gliAcyt_gliAnuc',1)

#where do I represent the conversion of Gli into GliA and GliR?(Parameter section?)
Rule('PTCH_binds_HH',PTCH(hh=None) + HH(ptch=None) | PTCH(hh=1) % HH(ptch=1), kf_ptch_hh, kr_ptch_hh)
Rule( 'Gli_binds_SUFU', Gli(sufu=None,state='U') + SUFU(gli=None) | Gli(sufu=1, state='U') % SUFU(gli=1), kf_gli_sufu, kr_gli_sufu)
Rule('Gli_phos_CK1', Gli(sufu=1, state='U') % SUFU(gli=1) + CK1() + PTCH(hh=None) >>
     Gli(sufu=1, state='P') % SUFU(gli=1) + CK1() + PTCH(hh=None), k_gli_phos_ck1)
Rule('Gli_phos_PKA', Gli(sufu=1, state='P') % SUFU(gli=1) + PKA() +PTCH(hh=None) >>
     Gli(sufu=1, state='PP') % SUFU(gli=1) + PKA() +PTCH(hh=None), k_gli_phos_pka)
Rule('Gli_phos_GSK3B', Gli(sufu=1, state='PP') % SUFU(gli=1) + GSK3B() + PTCH(hh=None) >>
     Gli(sufu=1, state='PPP') % SUFU(gli=1) + GSK3B() +PTCH(hh=None), k_gli_phos_gsk3b)
Rule('Gli_PPP_sufu', Gli(sufu=1, state='PPP') % SUFU(gli=1) >> Gli(sufu=None, state='PPP') + SUFU(gli=None), k_gli_ppp_sufu)
Rule('Gli_to_GliR', Gli(sufu=None, state='PPP') >> GliR(b=None, loc='cyt'),k_gli_gliR)
Rule('Gli_to_GliA_P', Gli(sufu=None, state='PPP') >> GliA(b=None, state='P', loc='cyt'), k_gli_gliA_P)
Rule('GliAp_degradation',GliA(b=None, state='P') >> None, k_gliA_deg)
Rule('Gli_GliA_U', Gli(sufu=None,state='U') >> GliA(b=None, state='U', loc='cyt'), k_gli_gliA_U)
Rule('GliR_Nu',GliR(b=None, loc='cyt')>> GliR(b=None, loc='nuc'), k_gliRcyt_gliRnuc)
Rule('GliA_Nu',GliA(b=None,state='U',loc='cyt') >> GliA(b=None, state='U', loc='nuc'), k_gliAcyt_gliAnuc)

#Rule('SMO_to_PC',SMO(pc=1,))

#Parameters for initials
Parameter('PTCH_init', 100)
Parameter('HH_init', 0)
Parameter('Gli_init', 100)
Parameter('SUFU_init', 100)
Parameter('CK1_init', 100)
Parameter('PKA_init', 100)
Parameter('GSK3B_init', 100)

#Initials
Initial(PTCH(hh=None), PTCH_init)
Initial(HH(ptch=None), HH_init)
Initial(Gli(sufu=None, state='U'), Gli_init)
Initial(SUFU(gli=None), SUFU_init)
Initial(CK1(), CK1_init)
Initial(PKA(), PKA_init)
Initial(GSK3B(), GSK3B_init)


Observable('GliR_nuc', GliR(loc='nuc'))
Observable('GliA_nuc', GliA(loc='nuc'))

# Observable('A_phos', A(state='P'))
# # Observable('A_bound, A(b=1') % B(a=1))
# Observable('A_bound', A(b=ANY))

# Presence of Hh Signaling Pathway

tspan=np.linspace(0,100,101)
sim=ScipyOdeSimulator(model,tspan,verbose=True)
result=sim.run()

plt.plot(tspan,result.observables['GliA_nuc'], lw=2,label='GliA_nuc')
plt.plot(tspan,result.observables['GliR_nuc'],lw=2,label='GliR_nuc')
plt.xlabel('time')
plt.ylabel('concentration')
plt.legend(loc='best')

plt.show()