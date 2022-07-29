from Dataloader import load_data
from Model import UNet
import torch
import torch.nn as nn
import matplotlib.pyplot as plt
import numpy as np
from random import choices
from geomloss import SamplesLoss
import time

def draw_samples(x, n, dtype=torch.FloatTensor):
    A = x
    xg, yg = np.meshgrid(
        np.linspace(0, 1, A.shape[0]),
        np.linspace(0, 1, A.shape[1]),
        indexing="xy",
    )

    grid = list(zip(xg.ravel(), yg.ravel()))
    dens = A.ravel() / A.sum()
    dots = np.array(choices(grid, dens, k=n))
    dots += (0.5 / A.shape[0]) * np.random.standard_normal(dots.shape)

    return torch.from_numpy(dots).type(dtype)

path = "/Users/martinkristiansen/Desktop/Simula_2022/Data/"
input = np.load(path+"input_images.npy")
target = np.load(path+"target_images.npy")

x = draw_samples(input[0], 100)

model = UNet()

loss_fn = SamplesLoss("sinkhorn", blur=0.01, scaling=0.9)
# x = model(x.unsqueeze(0))
# x = x[0, 0, :, :].detach()
x = x.clone().detach().requires_grad_(True)

def PointsInCircum(r,n=100):
    return [(0.4+np.cos(2*np.pi/n*x)*r,0.4+np.sin(2*np.pi/n*x)*r) for x in range(0,n+1)]

y = PointsInCircum(0.3, n=100)
y = torch.from_numpy(np.array(y)).type(torch.FloatTensor)

def display_samples(ax, x):
    x_ = x.detach().cpu().numpy()
    ax.scatter(x_[:, 0], x_[:, 1], s = 5)#, edgecolors="none")


# Euler Scheme
lr = 0.1
Nsteps = int(5 / lr) + 1
display_its = [int(t / lr) for t in [0, 0.25, 0.50, 1.0, 2.0, 5.0]]

t_0 = time.time()
plt.figure(figsize=(12, 8))
k = 1
for i in range(Nsteps):
    if i == 0:
        lr = 0
    else:
        lr = 0.1
    L = loss_fn(x, y)
    [g] = torch.autograd.grad(L, x)
    x.data -= lr *len(x)* g
    if i in display_its:  # display
        ax = plt.subplot(2, 3, k)
        k = k + 1
        plt.set_cmap("hsv")
        plt.scatter(
            [10], [10]
        )  # shameless hack to prevent a slight change of axis...

        display_samples(ax, y)
        display_samples(ax, x)

        ax.set_title("t = {:1.2f}".format(lr * i))

        plt.axis([0, 1, 0, 1])
        plt.gca().set_aspect("equal", adjustable="box")
        plt.xticks([], [])
        plt.yticks([], [])
        plt.tight_layout()

plt.show()


























#
#
# loss_values = []
# for epoch in range(num_epochs):
#     # zero the parameter gradients
#     optimizer.zero_grad()
#     # forward + backward + optimize
#     pred = model(x.unsqueeze(0))
#     loss = loss_fn(pred.squeeze(1).squeeze(0), y)
#     loss_values.append(loss.item())
#     loss.backward()
#     optimizer.step()
#     print(loss.item())
#
#     if epoch == 0:
#         continue
#     if loss_values[epoch] < loss_values[epoch-1]:
#         torch.save({
#                 'epoch': epoch,
#                 'model_state_dict': model.state_dict(),
#                 'optimizer_state_dict': optimizer.state_dict(),
#                 'loss': loss.item(),
#                 }, "Wasserstein.pt")
#
#     if epoch % 10 == 0:
#         # show prediction every 10th epoch
#         plt.cla()
#         plt.imshow(pred[0,0, :, :].detach())
#         plt.draw()
#         plt.pause(0.01)
