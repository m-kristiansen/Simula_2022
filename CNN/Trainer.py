import torch
import torch.nn as nn
import numpy as np
from Dataloader import load_data
from Model import UNet
import matplotlib.pyplot as plt

trainloader, testloader, train_data, test_data = load_data(split = 0.33, batch_size = 10, device = 'cuda:0')

def train(device, num_epochs, checkpoint):

    num_epochs = num_epochs
    model = UNet()

    #use this when training on the cluster
    device= torch.device(device)
    model = model.to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

    if checkpoint == True:
        checkpoint = torch.load("checkpoint.pt")
        model.load_state_dict(checkpoint['model_state_dict'])
        optimizer.load_state_dict(checkpoint['optimizer_state_dict'])

    model.train()

    loss_fn = nn.MSELoss()

    totalBatchLoss = np.zeros(num_epochs)
    for epoch in range(num_epochs):
        loss_values = []
        for i, (x, y) in enumerate(trainloader):
            # zero the parameter gradients
            optimizer.zero_grad()
            # forward + backward + optimize
            output = model(x).to(device)
            loss = loss_fn(output.squeeze(1), y)
            loss_values.append(loss.item())
            loss.backward()
            optimizer.step()



        totalBatchLoss[epoch] = sum(loss_values)
        print("loss at epoch {}/{}: {}".format(epoch, num_epochs, totalBatchLoss[epoch]))

        if epoch == 0:
            continue

        if totalBatchLoss[epoch] < totalBatchLoss[epoch-1]:
            torch.save({
                    'epoch': epoch,
                    'model_state_dict': model.state_dict(),
                    'optimizer_state_dict': optimizer.state_dict(),
                    'totalBatchLoss': totalBatchLoss,
                    }, 'checkpoint.pt')

train('cuda:0', 1000, True)
