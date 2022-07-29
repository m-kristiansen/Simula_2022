import torch
from Dataloader import load_data
from Model import UNet
import matplotlib.pyplot as plt
import numpy as np

trainloader, testloader, train_data, test_data = load_data(split = 0.33, batch_size = 64, device = 'cpu')

checkpoint = torch.load("checkpoint_cuda-0.pt", map_location = torch.device('cpu'))
losses = torch.load("losses_cuda-0.pt", map_location = torch.device('cpu'))

model = UNet()
model.load_state_dict(checkpoint['model_state_dict'])
trainloss = losses["TrainLoss"]
testloss = losses["TestLoss"]
best_epoch = checkpoint["best_epoch"]

model.eval()


def prediction_samples(test_data):
    inp, targ = test_data[0]
    pred = model(inp.unsqueeze(0))
    fig, axs = plt.subplots(2, 2)
    axs[0,0].imshow(pred[0, 0, :, :].detach())
    axs[0,0].set_title("Prediction")
    axs[0, 1].imshow(targ.detach())
    axs[0, 1].set_title("Target")
    inp, targ = test_data[1]
    pred = model(inp.unsqueeze(0))
    axs[1,0].imshow(pred[0, 0, :, :].detach())
    axs[1,1].imshow(targ.detach())
    plt.show()

prediction_samples(test_data)


plt.plot(range(len(trainloss[1:])), trainloss[1:], '.', label = 'Train')
plt.plot(range(len(testloss[1:])), testloss[1:], '.', label =  'Test')
plt.axvline(x = best_epoch , color = 'r', label = 'best_epoch')
plt.legend()
plt.show()
