"""
This file will contain the main method for computing our predictions of the image values with the graph as input
The input is given as a 128x128 binary matrix
the output is given as a 128x128 one channel image

The original UNet() uses channels = (1024, 512, 256, 128, 64) in both the encoder and decoder. This yields 31 million learnable parameters.
Since our input shape is (128, 128) and to reduce the number of parameters we restrict ourselfes to (256, 128, 64).
This yields 1.8 mill learnable parameters and an output shape before interpolation of (88, 88), meaning we retain most of our spatial size.
"""
import torch
import torch.nn as nn
import torch.nn.functional as F
import torchvision
import matplotlib.pyplot as plt
from torch.utils.data import Dataset, DataLoader

class Block(nn.Module):
    """ The block operation uses two 3x3 convolutions, each followed by a ReLU activation """
    def __init__(self, in_ch, out_ch):
        super().__init__()
        self.conv1 = nn.Conv2d(in_ch, out_ch, 3)#), padding = 'same')
        self.relu  = nn.ReLU()
        torch.nn.BatchNorm2d(out_ch)
        self.conv2 = nn.Conv2d(out_ch, out_ch, 3)#, padding = 'same')
        torch.nn.BatchNorm2d(out_ch)

    def forward(self, x):
        return self.relu(self.conv2(self.relu(self.conv1(x))))

class Encoder(nn.Module):
    """
    1. input size = (batch_size,1, 128, 128)
    block -> maxpool
    2. torch.Size([batch_size, 64, 62, 62])
    block -> maxpool
    3. torch.Size([batch_size, 128, 29, 29])
    block -> maxpool
    4. torch.Size([1, 256, 12, 12])
    """
    def __init__(self, chs=(1, 64, 128, 256)):
        super().__init__()
        self.enc_blocks = nn.ModuleList([Block(chs[i], chs[i+1]) for i in range(len(chs)-1)])
        self.pool       = nn.MaxPool2d(2)

    def forward(self, x):
        ftrs = []
        for i, block in enumerate(self.enc_blocks):
            x = block(x)
            ftrs.append(x)
            x = self.pool(x)
            #print("Encoder layer {}: {}".format(i, x.shape))
        return ftrs

class Decoder(nn.Module):
    """
    To concatenate between encoding and decoding we need the same dimensionality, the crop fixes this by
    centercropping the picture in the decoder with dimensions of the picture in the encoder
    In the decoder, block reduces the channel halves the channel dimensions.

    1. torch.Size([batch_size, 256, 12, 12])
    convTranspose -> crop -> concat -> block
    2. torch.Size([batch_size, 128, 46, 46])
    convTranspose -> crop -> concat -> block
    3. torch.Size([batch_size, 64, 88, 88])
    """

    def __init__(self, chs=(256, 128, 64)):
        super().__init__()
        self.chs        = chs
        self.upconvs    = nn.ModuleList([nn.ConvTranspose2d(chs[i], chs[i+1], kernel_size=2, stride=2) for i in range(len(chs)-1)])
        self.dec_blocks = nn.ModuleList([Block(chs[i], chs[i+1]) for i in range(len(chs)-1)])

    def forward(self, x, encoder_features):
        for i in range(len(self.chs)-1):
            x        = self.upconvs[i](x)
            enc_ftrs = self.crop(encoder_features[i], x)
            x        = torch.cat([x, enc_ftrs], dim=1)
            x        = self.dec_blocks[i](x)
            #print("Decoder layer {}: {}".format(i, x.shape))
        return x

    def crop(self, enc_ftrs, x):
        _, _, H, W = x.shape
        enc_ftrs   = torchvision.transforms.CenterCrop([H, W])(enc_ftrs)
        return enc_ftrs


class UNet(nn.Module):
    def __init__(self, enc_chs=(1, 64, 128, 256), dec_chs=(256, 128, 64), num_class=1):
        super().__init__()
        self.encoder     = Encoder(enc_chs)
        self.decoder     = Decoder(dec_chs)
        self.head        = nn.Conv2d(dec_chs[-1], num_class, 1)

    def forward(self, x):
        x = torch.unsqueeze(x, dim = 1)
        enc_ftrs = self.encoder(x)
        out      = self.decoder(enc_ftrs[::-1][0], enc_ftrs[::-1][1:])
        out      = self.head(out)
        out = F.interpolate(out, (128,128))
        return out


if __name__ == "__main__":
    #testing model, take 1 graph as input print the shape of the prediction and plot the prediction.
    #Note, this test is on an untrained model so the plot only gives you an idea of the impact of the model-
    #parameters on the graph.
    import numpy as np

    model = UNet()

    path = "/Users/martinkristiansen/Desktop/Simula_2022/Data/"
    input = np.load(path+"input_images.npy")

    x = torch.tensor(input[0, :, :]).unsqueeze(0).float() #convert
    print(x.shape)
    pred = model(x)
    print("output shape: ", pred.shape)
    plt.imshow(pred[0,0,:, :].detach())
    plt.show()
    pytorch_total_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print("total number of learnable parameters: ",pytorch_total_params)
