""" This file implements a simple euler scheme to check the perfomance of the loss function and backpropogation,
it also trains the model on one single image, this was used to choose the depth and the kernel sizes of the UNet.
"""
from Dataloader import load_data
from Model import UNet
import torch
import torch.nn as nn
import matplotlib.pyplot as plt
import numpy as np
path = "/Users/martinkristiansen/Desktop/Simula_2022/Data/"


def euler_scheme(x, y, num_epochs, loss_fn, lr, model):
    x = model(x.unsqueeze(0)).squeeze(0).squeeze(0) #make one prediction
    x = x.clone().detach().requires_grad_(True) # clear gradient history

    # Euler Scheme
    for epoch in range(num_epochs):
        L = loss_fn(x,y)
        print(L)
        L.backward()
        x.data -= lr * len(x)*x.grad
        print(L.item())
        plt.cla()
        plt.imshow(x.squeeze(0).detach())
        plt.pause(0.01)


def train_one_image(x, y, num_epochs, loss_fn, model, optimizer, use_checkpoint):
    checkpoint = torch.load("one_image.pt")

    if use_checkpoint == True:
        model.load_state_dict(checkpoint['model_state_dict'])
        optimizer.load_state_dict(checkpoint['optimizer_state_dict'])

    num_epochs = 100
    loss_values = []
    for epoch in range(num_epochs):
        # zero the parameter gradients
        optimizer.zero_grad()
        # forward + backward + optimize
        pred = model(x.unsqueeze(0))
        loss = loss_fn(pred[0, :, :], y.unsqueeze(0))
        loss_values.append(loss.item())
        loss.backward()
        optimizer.step()
        print("loss at epoch {}: {}".format(epoch, loss.item()))
        for name, param in model.named_parameters():
            # print norm of gradients in each layer
            print(name, param.grad.norm())

        if epoch == 0:
            continue
        if loss_values[epoch] < loss_values[epoch-1]:
            torch.save({
                    'epoch': epoch,
                    'model_state_dict': model.state_dict(),
                    'optimizer_state_dict': optimizer.state_dict(),
                    'loss': loss.item(),
                    }, "one_image.pt")

        if epoch % 10 == 0:
            # show prediction every 10th epoch
            plt.cla()
            plt.imshow(pred[0,0, :, :].detach())
            plt.draw()
            plt.pause(0.01)
    plt.show()


input = np.load(path+"input_images.npy")
target = np.load(path+"target_images.npy")

x = torch.from_numpy(input[0].astype(np.float32))
y = torch.from_numpy(target[0].astype(np.float32))

model = UNet()
optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
loss = nn.MSELoss()

#euler_scheme(x, y, 100, loss, 0.1, model)

train_one_image(x, y, 100, loss, model, optimizer, use_checkpoint = False)
