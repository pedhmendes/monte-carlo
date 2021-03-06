# Monte Carlo

Projects of Monte Carlo class.

All comments are in portuguese.

### Group Project - Ehrenfest Model
Done with my colleagues, Gabriela Oliveira and João Mesquita. The Jupyter Notebook can be found in this repo, the Colab link can be found [here](https://colab.research.google.com/drive/1wgrEgNMEajY78lugJqVN2RftX6Z6jVJh?usp=sharing).


### Ising Model
In Ising's folder can be found different forms to simulate the Ising Model. There's four codes: (1) simple, (2) without matrix, (3) with matrix and (4) bi array. The (1) simple code is the best one, is faster and can be easylly understand. I used the simple in my research, more info [here](https://github.com/pedhmendes/spin-systems-ic). The one (2) without matrix evaluate the energy variation in every step, is not the best way. (3) with matrix code I stored the topology of my system. The last one (4), bi array, uses bidimensional array, this one runs very slow. 

The ones with GSL flag in the end uses this C library, you can read more [here](https://github.com/pedhmendes/gsl).

There is also available the script to run *n* times, the execution times files, for the ones that did not used the gsl, and also an awk script to evaluate the mean and variance.

### Group Project - Potts Model
In Potts's folder can be found 4 codes for different simulation of Potts Model. There are two algorithms, metropolis and heat bath, and they are indicated in the code name. Also there are two versions, the ones with *if* and the ones with the *matrix*. The first ones are less optimized, they use a lot of if/else lines. The second ones are best, they do not have that much of if because they use a identity matrix that handle the comparissions we want to do. 

We did a wiki for this project, Caetano Slaviero and João Mesquita and I. You can find the wiki (in portuguese) in this [link](https://fiscomp.if.ufrgs.br/index.php/Modelo_de_Potts_2D). 

There is also available the script to run *n* times, for different *q* states and temperatures. And our usual library with MC functions can be found in the folder as well.

I would like to say thanks to Luis Latoski, graduate student, that help to develop and debug the codes.
