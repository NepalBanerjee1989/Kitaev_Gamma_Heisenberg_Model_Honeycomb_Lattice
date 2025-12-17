# Kitaev_Gamma_Heisenberg_Model_Honeycomb_Lattice
----------------------------------------------------------------------------
----------------------------------------------------------------------------
----------------------------------------------------------------------------
Author:-Nepal Banerjee,Department of Physics,University of Seoul,Seoul,South-Korea

Contact: Email Id :nepalbanerjee36@gmail.com/nb.uos1989@gmail.com\
         Mob Number:+91-9816075495/+91-9474172495/+91-8250773046
         

Copyright (c) 2025-2026 Nepal Banerjee. All rights reserved.

--------------------------------------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------
1)Here we have simulated the Kitaev-Gamma-Heisenberg Model for simulating the 2D  VdW material Alpha-RuCl3 using classical Monte-Carlo.Here we have predicted transition temperature T=6.5K, which is very close to experimental results reported on this material, which is 7-8K. Here, we also noticed the Classical Z2 spin liquid in the pure Kitaev limit.

2)Here, we use the parameters J1=0.1(FM), J2=0, J3=-0.2, K=1.0(FM), and Gamma=0.4(FM), which predict a  zigzag AFM ground state for this material. We confirmed the zigzag-AFM-like magnetic order with real-space spin configuration and with static spin-structure factor calculation, which confirmed that phase. Here we have seen the interplay of order by disorder mechanism, Heisenberg interaction,Kitaev interaction and gamma interaction and observed how this zigzag-AFM order developed in our system after tracking the microscopic configuration at low temperature. 

3)Here we have used the Fortran 77 language and the gfortran compiler for developing the project, and used Python for analysing the data.


We also simulated the Classical Z2 spin liquid using our simulation in the pure Kitaev limit,where K=4.0,J=Gamma=0.We notice a micro-vortex state from real space spin configuration and confirmed that an exotic phase of matter using the static spin structure factor calculation, which clearly revealed the existence of the classical Z2 spin liquid at temperature T=0.01.

4)We also observed classical frustrated induced U(1) gauge like spin liquid in purely Heisenberg limit.We use the exchange parameter J2=J3=J1/2 combination for simulating this frustrated spin liquid.Here we observed the spin liquid ground state with real space spin configuration and with a pinch point singularity atstatic structure factor S(q).

5)Here we also proposed a Hamiltonian which can coupled this two extreme limiting phase with a interpolating Hamiltonian with a tunning parameter "g" which can capture the phase transition from classical U(1) gauge spin liquid to classical Z2 flux spin liquid.

---------------------------------------------------------------------------
---------------------------------------------------------------------------
 Hamiltonian of our combined system:-
----------------------------------------------------------------------------
 H(g)=(1-g)*H_{Heisenberg} +g*H_{Kitaev}

---------------------------------------------------------------------------

Here "g" is a tunning parameter.When g=0 it is purely Heisenberg.When g=1 itis purely Kitaev like.In intermidate coupling g=1/2 we notice a Zigzag-AFM ground state.Our observation find good agreement with the following literature.

----------------------------------------------------------------------------
Reference:
----------------------------------------------------------------------------
1)Kitaev-Heisenberg Model on a Honeycomb Lattice:Possible Exotic Phases in Iridium Oxides A2IrO3.\
Authors:-Jiri Chaloupka,George Jackeli,and Giniyat Khaliullin

2)Quantum phase transition in Heisenberg-Kitaev model.\
Authors:-Robert Schaffer, Subhro Bhattacharjee, and Yong Baek Kim 

3)Spin-S Kitaev model: Classical ground states, order from disorder,
and exact correlation functions.\
Authors:-G.Baskaran, Diptiman Sen,and R. Shankar

4)Neutron scattering in the proximate quantum spin liquid a-RuCl3.\
Authors:-Arnab Banerjee,Jiaqiang Yan,Johannes Knolle,Craig A. Bridges,Matthew B. Stone,Mark D. Lumsden,David G. Mandrus,David A. Tennant,Roderich Moessner,Stephen E. Nagler

5)Monoclinic crystal structure of α-RuCl 3 and the zigzag antiferromagnetic ground state.\
Authors:-R.D.Johnson,S.C.Williams,A.A.Haghighirad,J.Singleton,V.Zapf,P.Manuel,I.I.Mazin,Y.Li,H.O.Jeschke,R.Valentı́,and R.Coldea

6)Classical Spin Liquid on the Maximally Frustrated Honeycomb Lattice.\
Authors:-J.Rehn,Arnab Sen,Kedar Damle,and R. Moessner

----------------------------------------------------------------------------


Future Plan:-
---------------------------------------------------------------------------
1)We will introduce further disorder effect like bond random disorder  and site impurity effect like quench disorder in this model.
2)We will apply our understanding in other model system.

------------------------------------------------------------------------------------------------------------------------------------------------------ 
Copyright (c) 2025-2026 Nepal Banerjee. All rights reserved.


--------------------------------------------------------------------------------------------------------------------------------------------------------

 Copy Right statement:                                    
--------------------------------------------------------------------------------------------------------------------------------------------------------
A copyright disclaimer is a notice stating that this content (code,text, images,video) is my intellectual property, preventing unauthorized use and clarifying ownership. it acts as a "No Trespassing" sign for my work,discouraging theft and helping with legal claims.For using others' work (like in videos), a "Fair Use" disclaimer citing Section 107 of theCopyright Act can justify limited use for education,commentary,etc.,but always credit the original owner.Github preserve and apply this copy right rules of intellectual properties strictly.With out permission of owner, using all this code,change and re-distribution of code is not allowed and consider as illegal,crime and punishable offence under the strict govt rules and regulation.
--------------------------------------------------------------------------------------------------------------------------------------------------------
