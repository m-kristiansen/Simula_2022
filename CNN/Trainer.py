import torch
import torch.nn as nn
import numpy as np
from Dataloader import load_data
from Model import UNet
import matplotlib.pyplot as plt


def train(device, num_epochs, checkpoint, testloader, trainloader):

    num_epochs = num_epochs
    model = UNet()

    #use this when training on the cluster
    device= torch.device(device)
    model = model.to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

    if checkpoint == True:
        checkpoint = torch.load("checkpoint_{}.pt".format(device), map_location = torch.device(device))
        model.load_state_dict(checkpoint['model_state_dict'])
        optimizer.load_state_dict(checkpoint['optimizer_state_dict'])


    loss_fn = nn.MSELoss()

    TrainLoss = []
    TestLoss = []
    for epoch in range(num_epochs):
        model.train()
        train_batch_losses = []
        for i, (x, y) in enumerate(trainloader):
            # zero the parameter gradients
            optimizer.zero_grad()
            # forward + backward + optimize
            output = model(x).to(device)
            loss = loss_fn(output, y.unsqueeze(1))
            train_batch_losses.append(loss.item())
            loss.backward()
            optimizer.step()

        with torch.no_grad():
            model.eval()
            test_batch_losses = []
            for i, (x_test, y_test) in enumerate(testloader):
                output = model(x_test).to(device)
                test_loss = loss_fn(output.squeeze(1), y_test)
                test_batch_losses.append(test_loss.item())

            TestLoss.append(np.mean(test_batch_losses))


        TrainLoss.append(np.mean(train_batch_losses))
        print("testloss: {} | trainloss: {} | at epoch {}/{}".format(TestLoss[-1], TrainLoss[-1], epoch, num_epochs))

        torch.save({
                'TrainLoss': TrainLoss,
                'TestLoss':TestLoss
                }, 'losses_{}.pt'.format(device))

        if epoch == 0:
            continue

        if TestLoss[-1] < np.min(TestLoss[0:-1]):
            torch.save({
                    'best_epoch': epoch,
                    'model_state_dict': model.state_dict(),
                    'optimizer_state_dict': optimizer.state_dict(),
                    }, 'checkpoint_{}.pt'.format(device))


if __name__ == "__main__":
    device = 'cuda:0'
    num_epochs = 100
    use_checkpoint = False
    batch_size = 10

    trainloader, testloader, train_data, test_data = load_data(split = 0.33, batch_size = batch_size, device = device)
    train(device, num_epochs, use_checkpoint, testloader, trainloader)
