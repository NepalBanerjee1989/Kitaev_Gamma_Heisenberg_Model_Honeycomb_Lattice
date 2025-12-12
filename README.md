# Kitaev_Gamma_Heisenberg_Model_Honeycomb_Lattice
--------------------------------------------------
Here we have simulated the Kitaev-Gamma-Heisenberg Model for simulating the 2D  VdW material Alpha-RuCl3 using classical Monte-Carlo. Here we have predicted transition temperature T=6.5K, which is very close to experimental results reported on this material, which is 7-8K. Here, we also noticed the Classical Z2 spin liquid in the pure Kitaev limit.
Here, we use the parameters J1=0.1(FM), J2=0, J3=-0.2, K=1.0(FM), and Gamma=0.4(FM), which predict a  zigzag AFM ground state for this material. We confirmed the zigzag-AFM-like magnetic order with real-space spin configuration and with static spin-structure factor calculation, which confirmed that phase. Here we have seen the interplay of order by disorder mechanism, Heisenberg interaction,Kitaev interaction and gamma interaction and observed how this zigzag-AFM order developed in our system after tracking the microscopic configuration at low temperature. 
Here we have used the Fortran 77 language and the gfortran compiler for developing the project, and used Python for analysing the data.
We also simulated the Classical Z2 spin liquid using our simulation in the pure Kitaev limit, where K=4.0, J=Gamma=0.
We notice a micro-vortex state from real space spin configuration and confirmed that an exotic phase of matter using the static spin structure factor calculation, which clearly revealed the existence of the classical Z2 spin liquid at temperature T=0.01.

We also observed classical frustrated induced U(1) gauge like spin liquid in purely Heisenberg limit.We use the J2=J3=J1/2 combination for simulating this frustrated spin liquid.Here we observed the spin liquid ground state and pinch point singularity at static structure factor S(q).

Here we also proposed a Hamiltonian which can coupled this two extreme limiting phase with a interpolating Hamiltonian with a tunning parameter "g" which can capture the phase transition from classical U(1) gauge spin liquid to classical Z2flux spin liquid.

Hamiltonian of our combined system:
-----------------------------------

H(g)=(1-g)*H_{Heisenberg} +g*H_{Kitaev}

When g=0 it is purely Heisenberg .When g=1 it is purely Kitaev like.


In intermidate coupling g=1/2 we notice a Zigzag-AFM ground state.Our observation find good agreement with the following literature :

Reference:-
----------
1)Kitaev-Heisenberg Model on a Honeycomb Lattice:Possible Exotic Phases in Iridium Oxides A2IrO3.
Author:-Jiri Chaloupka,George Jackeli,and Giniyat Khaliullin

2)Quantum phase transition in Heisenberg-Kitaev model.
Author:-Robert Schaffer, Subhro Bhattacharjee, and Yong Baek Kim 

3)Spin-S Kitaev model: Classical ground states, order from disorder,
and exact correlation functions.
Author:-G.Baskaran, Diptiman Sen,and R. Shankar

Future Plan:-
------------
1)We will introduce further disorder effect like bond random disorder  and site impurity effect like quench disorder in this model.
2)We will apply our understanding in other model system.
---------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------

Disclaimer Copy Right Rule:C,Nepal Banerjee,2025,"All Rights Reserved".
--------------------------

Disclaimer Copy Right Rules statement:-
-------------------------------------
A copyright disclaimer is a notice stating that your content (code,text, images,video) is your intellectual property, preventing unauthorized use and clarifying ownership. it acts as a "No Trespassing" sign for your work, discouraging theft and helping with legal claims. For using others' work (like in videos), a "Fair Use" disclaimer citing Section 107 of the Copyright Act can justify limited use for education, commentary, etc., but always credit the original owner.Github preserve and apply this copy right rules of intellectual properties strictly.With out permission of owner using the code and change and distribution of code is not allowed and consider as illegal,crime and punishable offence under the govt rules and regulation.

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
