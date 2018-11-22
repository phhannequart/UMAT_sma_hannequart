# UMAT_sma_hannequart
UMAT script for Abaqus, polycrystalline shape memory alloy model

For modelling a polycrystalline shape memory alloy wire
Authors : Philippe Hannequart (1,2), Michael Peigney(1), Jean-François Caron(1), Emmanuel Viglino(2)
(1) Universite Paris-Est, Laboratoire Navier (UMR 8205), CNRS, Ecole des Ponts ParisTech, IFSTTAR, 77455 Marne la vallee, France
(2) Arcora, Groupe Ingerop, Rueil-Malmaison, France
!
It refers to the model published in paper 
!
The number of state variables NSTATV in Abaqus must be twice the number of crystal orientations
NSTATV = 2 is a SMA monocrystal
!
The material properties in Abaqus must be set as follows :
PROPS(1) : Young's Modulus
PROPS(2) : Latent heat parameter
PROPS(3) : Reference temperature
PROPS(4) : Dissipation parameter for self-accomodated martensite 
PROPS(5) : Dissipation parameter for oriented martensite
PROPS(6) to PROPS(5 + NSTATV/2) are the values of the maximum phase transformation strains
