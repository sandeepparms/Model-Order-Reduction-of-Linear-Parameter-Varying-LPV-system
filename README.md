# Model Order Reduction of a Linear Parameter Varying (LPV) system

**In this project, we inspect the following spring-mass-damper system in which spring stifness is not constant but varies linearly according to a parameter**

<img src="https://user-images.githubusercontent.com/86435953/123509389-7fe27f80-d675-11eb-9652-e6d8a830bd65.PNG" width="700"  />

**Free body diagram of second mass with first and third mass attached to it**

![FBD](https://user-images.githubusercontent.com/86435953/123509433-c932cf00-d675-11eb-8062-f45ffd270338.PNG)

*Note: Files starting with name 'Multiple parameters..' are the main files and they have to be run. Other files are function files or evaluation files, install required toolbox and solvers before running the program*

In many engineering problems, a control engineer has to deal with a very large model with a huge number of states. This model is input for designing a controller. Two standard
control design techniques, H_inf and LQG result in a high order controller that is at least as complicated as the plant. This introduces undesired complexity in designing, analysis, and
synthesis of the controller. In almost in all practical applications, control engineer seeks a simpler solution. Model order reduction for linear parameter-varying (LPV) systems
is crucial in order to extend the applicability of the LPV framework to a wider range of applications. In this project work, a comparison of variants of the balanced truncation
method for LPV systems is carried out. The notion of an internally balanced realization is discussed. Detailed theory and computational procedure from Wood’s thesis (see literature) for arriving at the balanced truncation are discussed. To demonstrate the practicality of the discussed methods, a chained mass-spring-damper system is chosen and several LPV system formulations are made. The consequence of bounding the rate of parameter variation results in the complexity of the balanced system by introducing an exclusive dependency on the rate bound. Several combinations of mass-spring-damper system with different numbers of scheduling parameters are presented. The resulting balanced and truncated models are compared and discussed. The advantages and disadvantages of the LMI based approach are demonstrated.

## Time domain analyis of full order and truncated system

<img src="https://user-images.githubusercontent.com/86435953/123509221-70af0200-d674-11eb-828d-4d50682a43bc.png" width="500"  />

## Variation of singular values in LPV parameter space for 2 parameters

<img src="https://user-images.githubusercontent.com/86435953/123509225-76a4e300-d674-11eb-87d9-2d9f04293cd4.png" width="500"  />

## Bode plot: Full order system vs truncated system with 6 masses and 3 parametres

<img src="https://user-images.githubusercontent.com/86435953/123509253-9936fc00-d674-11eb-9635-bbaa9f56d0cb.png" width="500"  />
