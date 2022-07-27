import torch
from Dataloader import load_data
from Model import UNet
import matplotlib.pyplot as plt

trainloader, testloader, train_data, test_data = load_data(split = 0.33, batch_size = 10, device = 'cpu')

checkpoint = torch.load("checkpoint.pt", map_location = torch.device('cpu'))

model = UNet()
optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
model.load_state_dict(checkpoint['model_state_dict'])
optimizer.load_state_dict(checkpoint['optimizer_state_dict'])

model.eval()


inp, targ = test_data[0]

pred = model(inp.unsqueeze(0))
plt.imshow(pred[0, 0, :, :].detach())
plt.show()
plt.imshow(targ.detach())
plt.show()
