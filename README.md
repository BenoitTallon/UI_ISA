# UI_ISA intallation

Matlab app for the use of k_ISA_2D and k_ISA_3D functions.
UI_ISA can be intalled as standalone version (downloading matlab compiler) or as a Matalb App (App -> Install App).

# Running UI_ISA

In the drop down buttons : select materials available in library or enter custom properties selecting "custom".
Select the system geometry and click "start Calc."
Results are automatically displayed or can be save as a .mat file by clicking "Save results".

Frequency range and number of modes can be set in the top menu.

# k_ISA_2D and k_ISA_3D

[f,kISA,vph,vgr,ve,ltr,sigmaT] = Calc_k_ISA(Mat,Inc,N,fMAX,amoy,P,phi)

This function calculates the wavenumber in the frame of the ISA.

Inputs:

Mat = matrix
Inc = scatterers
N = number of modes (for the scattering function calculation)
fMAX = upper value for the frequency vector (MHz)
amoy = mean radius (mm)
P = polydispersity (%)
phi = filling fraction (%)

Outputs:

f = frequency vector (MHz)
kISA = wavenumber (mm-1)
vph = phase velocity (mm.μs^(-1))
vgr = group velocity (mm.μs^(-1))
ve = energy velocity (mm.μs^(-1))


    Note that if ‘Mat’ and ‘Inc’ inputs are characters; the program will take properties of the inclusions in the Choix_Materiaux.m function. 

For example: 
[f,kISA,cph,cgr,ve,ltr,sigmaT]=Calc_k_ISA('gel','verre',20,8,0.5,5,5);

    If you want to work with an unlisted material, you can replace the characters by a vector which contains the material properties [rho cL cT] (for fluids, take cT = 0).

For example, (glass beads in gel): 
[f,kISA,cph,cgr,ve,ltr,sigmaT]=Calc_k_ISA('gel', [2.5 5.60 3.4],20,8,0.5,5,5);

Or:
[f,kISA,cph,cgr,ve,ltr,sigmaT]=Calc_k_ISA([1 1.5 0], [2.5 5.60 3.4],20,8,0.5,5,5)


    Finally, the program does not work for solid matrix since the energy velocity is not computed in this case.



