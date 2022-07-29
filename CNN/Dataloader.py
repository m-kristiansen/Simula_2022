import torch
import numpy as np
from sklearn.model_selection import train_test_split
from torch.utils.data import Dataset, DataLoader
import torchvision.transforms.functional as TF
import sys

# Convert data to torch tensors
class Data(Dataset):
    def __init__(self, x, y, device):
        self.x = torch.from_numpy(x.astype(np.float32)).to(device)
        self.y = torch.from_numpy(y.astype(np.float32)).to(device)
        self.len = self.x.shape[0]

    def transforms(self, x, y):
        # Random horizontal flipping
        if torch.rand(1) > 0.5:
            x = TF.hflip(x)
            y = TF.hflip(y)

        # Random vertical flipping
        if torch.rand(1) > 0.5:
            x = TF.vflip(x)
            y = TF.vflip(y)
        return x, y


    def __getitem__(self, index):
        input_image, target_image = self.x[index], self.y[index]
        x, y = self.transforms(input_image, target_image)
        return x,y

    def __len__(self):
        return self.len

def load_data(split, batch_size, device):
    if device == 'cpu':
        path = '/Users/martinkristiansen/Desktop/Simula_2022/Data/'
        input = np.load(path+"input_images_10k.npy")
        target = np.load(path+"target_images_10k.npy")

    if device == 'cuda:0':
        input = np.load("input_images_10k.npy")
        target = np.load("target_images_10k.npy")

    #split
    x_train, x_test, y_train, y_test = train_test_split(input, target, test_size=split)

    # Use troch.utils functionality to initiate data
    train_data = Data(x_train, y_train, device)
    train_dataloader = DataLoader(dataset=train_data, batch_size=batch_size, shuffle=True)
    test_data = Data(x_test, y_test, device)
    test_dataloader = DataLoader(dataset=test_data, batch_size=batch_size, shuffle=True)

    return train_dataloader, test_dataloader, train_data, test_data

if __name__ == "__main__":
    #check if its working properly
    test_data = load_data(0.33, 10)[1]
    for batch, (x, y) in enumerate(test_data):
        print(f"Batch: {batch+1}")
        print(f"x shape: {x.shape}")
        print(f"y shape: {y.shape}")
