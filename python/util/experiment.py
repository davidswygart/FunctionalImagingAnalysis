# import pyximport
# pyximport.install(setup_args={"script_args" : ["--verbose"]})
# if __name__ == "__main__":
#     from . import preprocess
#     import exp
#     from utils import intervals, chirp_pattern
# else:
from . import preprocess
from . import exp
from . import register
from .utils import intervals, chirp_pattern

import numpy as np
import pandas as pd

from typing import List

from scipy.interpolate import interp1d

import os
import suite2p

from PIL import Image





#TODO: we should accept a list of raw files
#we can save the flips as an offset from the first file, maybe?
#instead of saving the frame index for every pixel, we can just save the index for each flip?

#TODO: separate out the common experiment code from the chirp/spot/bar code?

class Experiment:
    def __init__(self, timing_file: str, raw_files: List[str], bin_path: str, artifact_file: str, stack_file: str, mask_file: str, dark_file: str, func_channel = 1, anat_channel=0, planes=None):
        self.timing_file = timing_file
        self.raw_files = raw_files
        self.bin_path = bin_path
        self.func_channel = func_channel

        artifact = np.load(artifact_file).astype(np.double)
        artifact[artifact==0] = np.nan



        epochs = preprocess.flatten_epochs(preprocess.extract_parameters(self.timing_file))
        
        for param in ['LED','intensity']:
            epochs[param] = (epochs.trial_type == b'field') * epochs[f'spot{param[0].upper()}{param[1:]}'] + (epochs.trial_type == b'chirp') * epochs[f'chirp{param[0].upper()}{param[1:]}']  + (epochs.trial_type == b'bars') * epochs[f'bar{param[0].upper()}{param[1:]}']
        
        
        epochs['end_time'] = np.uint64(epochs.preTime + epochs.stimTime + epochs.tailTime) * 10000 + epochs.start_time
        #NOTE: pre/stim/post times are in units of milliseconds. Multiplying by 10k gives same delta as .NET
        epochs = epochs.sort_values('start_time')
        epochs['flips'] = None
        epochs['frame_indices'] = None
        
        mask = Image.open(mask_file)
        mask_shape = mask.size
        mask = np.asarray(mask.getdata()).reshape(mask_shape)

        #TODO: we need to increment flips for each file?
        
        for raw_file in self.raw_files:
            raw_file_root = os.path.splitext(os.path.split(raw_file)[-1])[0]
            print(raw_file_root)

            # 1) Read imaging data
            raw, props = preprocess.read_raw_file(raw_file, None)
            npix = props['frame_shape'][0] * props['frame_shape'][1]
            
            # 2) Align to stimulus/projector
            fi = exp.align_frames_to_epochs(props['time_stamps'], epochs.start_time.to_numpy(), epochs.end_time.to_numpy())
            trigger, flips = preprocess.get_trigger_times(raw, props)

            for i, (start, end) in enumerate(fi):
                if start<0 or end<0:
                    continue
                #TODO: we assume a projector frame rate of 60?
                tfi, tflips, *_= exp.cluster_trigger_events(trigger, flips, (start-16)*npix, (end+16)*npix, props['pixel_rate'], epochs.iat[i, epochs.columns.get_loc('frameRate')])
                epochs.iat[i, epochs.columns.get_loc('flips')] = tflips
                epochs.iat[i, epochs.columns.get_loc('frame_indices')] = tfi

            iepochs = epochs[np.logical_and(fi[:,0]>=0, fi[:,1]>=0)]
        
            spots, chirps, bars = exp.time_events(iepochs)


            # 3) Reject stimulus artifact
            ievents = pd.concat((chirps,bars,spots), ignore_index=True, join='inner').join(iepochs, on='epoch_id')
            ievents = ievents[ievents.frame_indices.apply(max) <= 2098]
            self.events = ievents
            ievents.to_parquet(os.path.join(bin_path, raw_file_root + '_events.parquet'))           

            if not os.path.exists(os.path.join(bin_path, raw_file_root + '_functional.npy')) or not os.path.exists(os.path.join(bin_path, raw_file_root + '_anatomy.bin')):
                
                self.anat = raw[:,anat_channel].reshape((-1,*props['frame_shape']))
                self.func = exp.mask_artifact(raw[:,func_channel], ievents, props, artifact).reshape((-1,*props['frame_shape']))
                
                # 4) Interpolate artifact in anatomy channel
                for i in range(props['frame_shape'][0]):
                    for j in range(props['frame_shape'][1]):
                        n = np.isnan(self.func[:,i,j])
                        self.anat[n,i,j] = interp1d(np.nonzero(~n)[0], self.anat[~n,i,j], fill_value='extrapolate')(np.nonzero(n)[0])

                # 5) Save epoched functional and anatomical data

                if not os.path.exists(self.bin_path):
                    os.mkdir(self.bin_path)

                np.save(os.path.join(bin_path, raw_file_root + '_functional.npy'), self.func)
                anat_file = suite2p.io.BinaryFile(Lx=props['frame_shape'][1], Ly=props['frame_shape'][0], n_frames = len(self.anat), filename=os.path.join(bin_path, raw_file_root + '_anatomy.bin'))
                anat_file[np.arange(len(self.anat))] = self.anat
            else:
                anat_file = suite2p.io.BinaryFile(Lx=props['frame_shape'][1], Ly=props['frame_shape'][0], filename=os.path.join(bin_path, raw_file_root + '_anatomy.bin'))
                self.func = np.load(os.path.join(bin_path, raw_file_root + '_functional.npy'))
                
            # 6) register
            if os.path.exists(os.path.join(bin_path, raw_file_root + '_anatomy_reg.npy')):
                print('Already registered')
            else:
                register.register(
                    anat_file,
                    stack_file,
                    props,
                    planes=planes,
                )
            del anat_file

            # 7) segment
            
            if len(spots) and not os.path.exists(os.path.join(bin_path, raw_file_root + '_spots.parquet')):
                si = spots.copy()
                si.index = si.index + si.epoch_id*28 #nspotsperepoch
                spots_df = preprocess.segment(os.path.join(bin_path, raw_file_root), props, si, mask, dark_file, func_channel = 1, data = self.func)
                spots_df.join(si[['cx','cy']], on='trial').to_parquet(os.path.join(bin_path, raw_file_root + '_spots.parquet'))
                del spots_df

            if len(chirps) and not os.path.exists(os.path.join(bin_path, raw_file_root + '_chirps.parquet')):
                chirps_df = preprocess.segment(os.path.join(bin_path, raw_file_root), props, chirps.set_index('epoch_id'), mask, dark_file, func_channel = 1, data = self.func)
                chirps_df.to_parquet(os.path.join(bin_path, raw_file_root + '_chirps.parquet'))
                del chirps_df

            if len(bars) and not os.path.exists(os.path.join(bin_path, raw_file_root + '_bars.parquet')):
                bi = bars.copy()
                bi.index = bi.index + bi.epoch_id*10 #nbarsperepoch
                bars_df = preprocess.segment(os.path.join(bin_path, raw_file_root), props, bi, mask, dark_file, func_channel = 1, data = self.func)
                bars_df.join(bi[['theta']], on='trial').to_parquet(os.path.join(bin_path, raw_file_root + '_bars.parquet'))
                del bars_df

        # 8) interpolate?
        

if __name__ == "__main__":
    exp_name = "071323B"

    # symphony_file = f"C:\\Users\\zfj\\data\\{exp_name}\\{exp_name}.h5"
    # # raw_file = r"C:\Users\zfj\data\021923B\n1337_p1241_00007.tif"
    # # raw_file = r"C:\Users\zfj\data\021923B\n1337_p1241_00008.tif"
    # raw_file_root = r"region1_00001"

    # raw_file = f"C:\\Users\\zfj\\data\\{exp_name}\\{raw_file_root}.tif"
    # bin_path =f"C:\\Users\\zfj\\data\\{exp_name}\\func"
    # mask_file = f"C:\\Users\\zfj\\data\\{exp_name}\\071323B_region1_stack_red_max_vitreousmasked_cp_masks.png"

    # stack = f"C:\\Users\\zfj\\data\\{exp_name}\\071323B_region1_stack_red_max_vitreousmasked.tif"

    # artifact_file = f"C:\\Users\\zfj\\code\\functionalimaginganalysis\\python\\artifact_lut.npy"
        

    # e = Experiment(symphony_file, [raw_file], None, artifact_file=artifact_file)

    # exp_name = "081023B"
    artifact_file = "C:\\Users\\zfj\\code\\functionalimaginganalysis\\python\\artifact_lut.npy"
    events_file = "C:\\Users\\zfj\\data\\081023B\\func\\region1_927nm_00001_events.parquet"

    artifact = np.load(artifact_file).astype(np.double)
    artifact[artifact==0] = np.nan

    ievents = pd.read_parquet(events_file)
    
    
    raw_file = "C:\\Users\\zfj\\data\\081023B\\region1_927nm_00001.tif"
    raw, props = preprocess.read_raw_file(raw_file, None)
            

            # e.events.iloc[26:27]
    #a particular bar trial...
    ma = exp.mask_artifact(raw[:,1], ievents.iloc[26:27], props, artifact)
    np.save("C:\\Users\\zfj\\data\\081023B\\func\\region1_927nm_00001_functional2.npy", ma)