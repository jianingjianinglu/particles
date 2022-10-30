# particles
The particles experience Lennard-Jones potential. When the distance between them is shorter, they repulse each other. 
When the distance between them is longer, they attract each other. 

Since the system is isolated, there’re no external particles acting on the system. 
The potential energy on every particle by another is V/2, where V is defined using the Lennard-Jones formula.
The kinetic energy on every particle is 1/2*m*v^2.

I calculated the ke and pe of one particle exterted by another particle, and then added up the ke and pe on every particle from every other particle in the system.

The total energy should be conserved.
Some peaks of deviation did appear. The max difference in the values is around 0.2 units of energy, so considering the peaks, the error is around 2%.
