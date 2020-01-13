# Arduino_Constrained_MPC_Library
This is a compact Constrained (linear) Model Predictive Control (MPC) library for Teensy4.0/Arduino system (or real time embedded system in general).
- It's not using Eigen (small source code - more simple to understand).
- It's not using C++ Standard Library/std (for embedded consideration).
- If you set `SYSTEM_IMPLEMENTATION` to `SYSTEM_IMPLEMENTATION_EMBEDDED_NO_PRINT` in `konfig.h`, the code is platform agnostic (not using any library beside these C header files: `stdlib.h`, `stdint.h`, and `math.h`).
- There's no malloc/new/free dynamic memory allocation for real time application (but using heavy stack local variables, so you need to run it through static memory analyzer if you are really concerned about implement this in mission critical hard real time application).

# The Background
This is the continuation of my previous repo ([Arduino_Unconstrained_MPC_Library](pronenewbits/Arduino_Unconstrained_MPC_Library#)) so you should read it before continue. As with my other repo, the main goal is for the student to learn the MPC concept (I've made decision to sacrifice speed to get best code readability I could get) while still capable of tackling real-time control system implementation (the code is computed in **3 ms to 7 ms**! See *Some Benchmark* section below).

To recap, the MPC formula derivation can be described as (I'm using Jan Maciejowski's *Predictive Control with Constraints* as reference, great book btw) :
![MPC derivation](Penurunan.png "Click to maximize if the image rescaling make you dizzy")

We have 3 type of constraints:
1. <img src="http://latex.codecogs.com/gif.latex?\Delta&space;U\left&space;(k&space;\right&space;)" border="0"/> constraint (i.e. <img src="http://latex.codecogs.com/gif.latex?\Delta&space;U\textsubscript{min}\leq&space;\Delta&space;U\left&space;(k&space;\right&space;)\leq&space;\Delta&space;U\textsubscript{max}" border="0"/>), or the actuator slew rate constraint.
2. <img src="http://latex.codecogs.com/gif.latex?u\left&space;(k&space;\right&space;)" border="0"/> constraint (i.e. <img src="http://latex.codecogs.com/gif.latex?u\textsubscript{min}\leq&space;u\left&space;(k&space;\right&space;)\leq&space;u\textsubscript{max}" border="0"/>), or the actuator output constraint.
3. <img src="http://latex.codecogs.com/gif.latex?z\left&space;(k&space;\right&space;)" border="0"/> constraint (i.e. <img src="http://latex.codecogs.com/gif.latex?z\textsubscript{min}\leq&space;z\left&space;(k&space;\right&space;)\leq&space;z\textsubscript{max}" border="0"/>), or the system output constraint.

**Note**: To make the explanation simple, in this implementation I use full rank constraints. But in the real implementation, you don't need to (e.g. for <img src="http://latex.codecogs.com/gif.latex?u\left (k \right )" border="0"/> you can only implement the minimum constraint (without the maximum constraint) <img src="http://latex.codecogs.com/gif.latex?u\textsubscript{min}\leq&space;u\left&space;(k&space;\right&space;)&space;;&space;u\left&space;(k&space;\right&space;)&space;\epsilon&space;R\textsuperscript{M}" border="0"/>). Or for example you have 3 input, you can only constraining the second input <img src="http://latex.codecogs.com/gif.latex?u\textsubscript{2}\left&space;(k&space;\right&space;)\leq&space;u\textsubscript{2,max}" border="0"/> and set free the other.

The constraints formulation can be described as:
![MPC constraints derivation](Constraints.png "Click to maximize if the image rescaling make you dizzy")
<p align="center"><small>welp, I guess that's all :3</small></p>

&nbsp;

Then we can describe the full formulation for Constrained MPC:
![MPC full formulation](Formulasi_Permasalahan_MPC.png "Click to maximize if the image rescaling make you dizzy")

*note:* You don't need to implement the full contraints matrix. Actually it is preferred to implement hard-contraints as little as possible to ensure maximum feasible search space.

# The Quadratic Programming Solver
In this implementation, I use (one of many) QP solver called [Active Set](https://en.wikipedia.org/wiki/Active-set_method). The big idea of active set algorithm is searching the optimal `x` value by solving the QP problem (with inequality constraints) as QP problem with equality constraints (EQP). The Active Set algorithm used in this implementation can be described as: 
![Active Set Algorithm](ActiveSet.png "Click to maximize if the image rescaling make you dizzy")

# Wrap it up
The constrained MPC then can be described as:
![Constrained MPC Algorithm](Kalkulasi.png "Click to maximize if the image rescaling make you dizzy")


# How to Use
Bla bla bla lorem ipsum bla bla


# Some Benchmark
Bla bla bla lorem ipsum bla bla



# Closing Remark
The matrix.h library's code documentation is still in Indonesian, but I plan to translate it into English soon (stay tuned!). In the meantime, it will be nice if you can test & validate my result or inform me if there are some bugs you encounter along the way! (or if you notice some grammar error in the documentation).

I published the code under CC0 license, effectively placed the code on public domain. But it will be great if you can tell me if you use the code, for what/why. That means a lot to me and give me motivation to expand the work (⌒▽⌒)
