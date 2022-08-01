# Simula 2022

This is a short summery of the summer project at the Simula Research Laboratory which consisted of mainly three parts

1. Random graph generation             -  Graphs/ 
2. Solving the coupled 2D/1D system   -  Solver/
3. Convolutional neural network         -  CNN/

The first part consist of producing random vascular graphs where the points/nodes are generated using a rapid random tree algorithm and connected such that it resembles the spread of real vasculature. 

![Graphs](https://github.com/m-kristiansen/Simula_2022/blob/main/Figures/4_graphs.png)

In the second part we want to study how different graphs describes different concentrations in the surrounding tissue. The following equations governs the coupled 2D/1D system 

$$
\begin{align*}
-\Delta u &= f   \hspace{1cm} &\text{in} \hspace{0.3cm}\Omega\\
[\[u]\] &= 0 \hspace{1cm} &\text{on}\hspace{0.3cm} \Gamma \\
-\Delta \hat{u} +  [\[-\nabla u \cdot \nu ]\] &= \hat{f} \hspace{1cm} &\text{on} \hspace{0.3cm}\Gamma \\
u - \hat{u} &= g \hspace{1cm} &\text{on}\hspace{0.3cm} \Gamma 
\end{align*}
$$

To solve this system we use the functionality of gmsh to generate our finite element mesh and the functionality of Gridap to enforce our boundary conditions and solving the weak formulation of our system on the mesh. The boundary conditions are zero Neumann on the bounding box walls, one Dirichlet on the root of the graph and 1/(1+x) Dirichlet on the end tips of the graph, where x is the distance from the root to the tip along the graph. All the source terms are set to zero. 

<img src="https://github.com/m-kristiansen/Simula_2022/blob/main/Figures/solution.png" width="500">

The last part consist of using ML to create a function taking the 2D binary matrix representation of the graph as input and prediciting an interpolation of the solution from 2., i.e. a single channel image. The chosen network is pure convolutional and since we want the input size = output size we use a UNet architecture. Below are four examples of inputs and targets for our model. These are also randomly flipped to make the model learn the more generalized features. 

Input to the network:
<img src = "https://github.com/m-kristiansen/Simula_2022/blob/main/Figures/4_binary_masks.png" width = "1000">

Target: 
<img src = "https://github.com/m-kristiansen/Simula_2022/blob/main/Figures/4_interpolations.png" width = "1000">

Our reference paper (https://www.sciencedirect.com/science/article/abs/pii/S0045793020302164) consist of a much easier machine learning problem, predicting the drag coefficient with 300,000 learnable parameteres used 12,000 input shapes. We want to predict 128x128 pixel values using 1,8 million learnable parameters, so clearly we need more data. However, I think the figures below gives an illustration of the learnability of the network. 

<p float="left">
<img src = "https://github.com/m-kristiansen/Simula_2022/blob/main/Figures/losses_1000.png" width = "400"/>
<img src = "https://github.com/m-kristiansen/Simula_2022/blob/main/Figures/pred_1000.png" width = "400"/>
</p>
<p float="left">
<img src = "https://github.com/m-kristiansen/Simula_2022/blob/main/Figures/losses_10k.png" width = "400"/>
<img src = "https://github.com/m-kristiansen/Simula_2022/blob/main/Figures/pred_10k.png" width = "400"/>
</p>

Sadly the project ends now, but I will list a few improvements that should be done for potential future endeavours. 

- More data, data augmentation can be applied on the 10k already produces different graphs. 
- Hyperparameter tuning
- Cross validation 
- Explore the ability to let the graph to connect with the box edges, and also allow the graphs to have roots along the box sides, not just the corners. 

There is also this thing called Wasserstein distance or earths mover distance which measures how far apart two sets of distribution's are. It would be fun to explore if this could be used as a loss function for our network. 
https://www.kernel-operations.io/geomloss/index.html has an pytorch compatible implementation measuring distances of points clouds. The question is if this also can be applied to one channel images. Here's a small illustration of the geomloss functionality:

<img src = "https://github.com/m-kristiansen/Simula_2022/blob/main/Figures/Wasserstein.png" width = "400"/>

Lastly I want to mention this awesome article https://www.nature.com/articles/s41598-021-85434-9 which inspire a way to produce much more realistic graphs.  

Main files:

Graphs/
- RRT.jl    - Generates RRT points and connects them. 
- graphs.jl - Creates bounding box for graph and generates mesh. Provides also different tags needed to solve the coupled system with spesific BC's.

Solver/
- Solver.jl - Solves the coupled system. 
- Interpolate.jl - Interpolates the solution on given partition. Used to generate the target images. 
- Binary_mask.jl - Creates binary mask of graph on given partition. Used to generate the input images. 
- get_data.jl - Combines the above to to generate data in iterative loop.

CNN/
- Dataloader.py - Uses Pytorch functionality to load data with given test/train split and batch size. 
- Model.py - Unet implementation. 
- Trainer.py - Train the network with given loss and optimizer. 
- Eval.py - Evaluate best model and test/train loss. 

