import suite2p
import re
import numpy as np
import argparse
from scipy.interpolate import interp1d
import os

def main(inpath, outpath, channel=0, load_in_memory=True, resample_uniform=True):

    path,file_start = os.path.split(inpath)
    sd = []
    for file in os.listdir(path):
        if not file.startswith(file_start):
            continue

        stack = suite2p.io.tiff.ScanImageTiffReader(os.path.join(path, file))
        md = stack.metadata()
            
        r = re.search("SI.hStackManager.framesPerSlice = (\d*)(Inf)?", md)
        if r is None: #TODO: for now, assumes a single frame/plane
            print(md)
            raise Exception("File is not a scanimage tiff")
        else:
            if not len(r.groups()[0]):
                n_frames_per_plane = int(re.search("SI.hScan2D.logAverageFactor = (\d+)", md).groups()[0])
            else:
                n_frames_per_plane = int(r.groups()[0]) // int(re.search("SI.hScan2D.logAverageFactor = (\d+)", md).groups()[0])

            stack_n_chans = len(re.search("SI.hChannels.channelSave = \[((?:\d+[ ;]?)+)\]", md).groups()[0].replace(' ',';').split(';'))
            
        stack_size = stack.shape()[-2:]
        if load_in_memory:
            stack_data = stack.data().reshape((-1,n_frames_per_plane,stack_n_chans,*stack_size))[:,:,channel].mean(axis=1)
        else:
            n_pages = n_frames_per_plane//stack_n_chans
            
            n_planes = stack.shape()[0]//n_pages

            stack_data = np.empty((n_planes,*stack_size), dtype=np.int32)

            for i in range(n_planes):
                stack_data[i] = stack.data(beg=n_pages*i, end=n_pages*(i+1)).reshape((n_frames_per_plane,stack_n_chans,*stack_size))[:,channel].mean(axis=0)

        po = suite2p.registration.bidiphase.compute(stack_data)
        if po != 0:
            suite2p.registration.bidiphase.shift(stack_data, po)

        if resample_uniform and (re.search("SI.hScan2D.uniformSampling = true", md) is not None):
            tff = float(re.search("SI.hScan2D.fillFractionTemporal = (0.\d+)", md).groups()[0])
            fov = [float(x) for x in re.split(r'[\s;]+', re.search("SI.hRoiManager.imagingFovUm = \[((?:\s*;*-?\d+.?\d*)*)\]", md).groups()[0])]
            
            stack_data = resample_uniform_scan(stack_data, fov, tff)
        
        sd.append(stack_data)

    suite2p.io.tiff.save_tiff(np.concatenate(sd, axis=0), outpath)

def resample_uniform_scan(data, fov, tff):
    yres = data.shape[1] / (fov[5] - fov[1]) #pixel per um
    xs = np.round((fov[2] - fov[0]) * yres).astype(int) # pixels

    res = np.empty((data.shape[0], data.shape[1], xs), dtype=np.int32)

    xi = np.sin(np.linspace(-np.pi/2 * tff, np.pi/2 * tff, data.shape[2])) # the sampled points
    qi = np.linspace(-tff, tff, xs) # the resampled points
    for z in range(data.shape[0]):
        for y in range(data.shape[1]):
            res[z,y] = interp1d(xi, data[z,y])(qi)
    return res

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog = 'ScanImage Stack processor',
        description = 'This program extracts a given channel from a ScanImage tiff, averages frames within a plane, and saves the result.',
    )
    parser.add_argument('input_file')
    parser.add_argument('output_file')
    parser.add_argument('-c','--channel', type=int, default=0)
    parser.add_argument('-f','--fast', type=bool, default=True)
    parser.add_argument('-u','--uniform', type=bool, default=True)

    args = parser.parse_args()
    main(args.input_file, args.output_file, channel = args.channel, load_in_memory = args.fast, resample_uniform=args.uniform)
