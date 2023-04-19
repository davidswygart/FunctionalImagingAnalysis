#general
import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Any
import os

#symphony
import h5py

#scanimage
import suite2p
import datetime
import re

if __name__ == "__main__":
    import register
else:
    from . import register

from skimage.transform import resize

#triggering
from sklearn.cluster import KMeans
from scipy.signal import find_peaks

from numba import jit, vectorize


#constant
DOTNET_EPOCH = datetime.datetime(1, 1, 1)
US = datetime.timedelta(microseconds = 1)

class Protocol():
    def __init__(self):
        self.epoch_blocks = []

class SpotsAndChirpAndBarsProtocol(Protocol):
    def __call__(self, name, g):
        if 'protocolID' in g.attrs and g.attrs['protocolID'] == b'sa_labs.protocols.stage.SpotFieldAndChirpAndBars': 
            attrs = g['protocolParameters'].attrs
            
            n_spots = int(attrs['frameRate'] * 35 /(attrs['spotPreFrames'] + attrs['spotStimFrames'] + attrs['spotTailFrames']))
            n_bars = 10

            emp = np.ones(len(g['epochs'])) * np.nan
            epochs = {
                'start_time': np.empty_like(emp, dtype=np.uint64),
                'trial_type': np.empty_like(emp, dtype=object),
                'cx': np.empty((len(g['epochs']), n_spots)) * np.nan,
                'cy': np.empty((len(g['epochs']), n_spots)) * np.nan,
                'theta': np.empty((len(g['epochs']), n_bars)) * np.nan,
            }   
            for i,epoch in enumerate(g['epochs']):
                epochs['start_time'][i]=g['epochs'][epoch].attrs['startTimeDotNetDateTimeOffsetTicks']
                epochs['trial_type'][i]=g['epochs'][epoch]['protocolParameters'].attrs['trialType']
                if epochs['trial_type'][i] == b'field':
                    epochs['cx'][i,:]=g['epochs'][epoch]['protocolParameters'].attrs['cx']
                    epochs['cy'][i,:]=g['epochs'][epoch]['protocolParameters'].attrs['cy']
                elif epochs['trial_type'][i] == b'bars':
                    epochs['theta'][i,:]=g['epochs'][epoch]['protocolParameters'].attrs['theta']
                elif epochs['trial_type'][i] == b'chirp':
                    pass
                
            self.epoch_blocks.append({**attrs, 'epochs':epochs})

def extract_parameters(timing_file: str) -> List[Dict[str, List]]:
    '''
    '''

    symphony_file = h5py.File(timing_file,"r")
    p = SpotsAndChirpAndBarsProtocol()
    symphony_file.visititems(p)

    return p.epoch_blocks

def flatten_epochs(epoch_blocks: List[Dict[str, Any]], block_params_to_keep: List[str] = None) -> pd.DataFrame:
    ''' Organize the metadata as a `structure of arrays`, with each index represeting one epoch
    By default, metadata pertaining to an entire block of epochs will be appended to each epoch
    Optionally specify which block metadata to keep
    '''
    if block_params_to_keep is None:
        block_params_to_keep = list(epoch_blocks[0].keys())
        block_params_to_keep.remove('epochs')

        dtype = [
            *[
                (epoch_blocks[0]['epochs'][k].dtype
                if type(epoch_blocks[0]['epochs'][k]) is np.ndarray
                and type(epoch_blocks[0]['epochs'][k][0]) is not np.ndarray
                else object)
                for k in epoch_blocks[0]['epochs'].keys()],
            *[epoch_blocks[0][k].dtype for k in block_params_to_keep]]
        keys = [*[k for k in epoch_blocks[0]['epochs'].keys()], *[k for k in block_params_to_keep]]
        dtype_d = {k:(str(t) if np.issubdtype(t,np.number) else object) for k,t in zip(keys,dtype)}


        return pd.DataFrame.from_dict({
                **{k: [v for eb in epoch_blocks for v in eb['epochs'][k]] for k in epoch_blocks[0]['epochs'].keys()},
                **{k: [eb[k] for eb in epoch_blocks for _    in eb['epochs']['trial_type']] for k in block_params_to_keep},
            }, orient='columns').astype(dtype_d)

def get_metadata(raw_data_file: str) -> str:
    return suite2p.io.tiff.ScanImageTiffReader(raw_data_file).metadata()

def read_raw_file(raw_data_file: str, raw_binary_path: str, functional_channel: int=1, anatomy_channel: int=0) -> Dict[str, Any]:
    ''' Read the scanimage file(s) for timestamps, sample rate
    '''

    props = {}     

    tif = suite2p.io.tiff.ScanImageTiffReader(raw_data_file)
    md = tif.metadata()

    props['frame_rate'] = float(re.search("SI.hRoiManager.scanFrameRate = (\d+.?\d*)", md).groups()[0])
    props['line_rate'] = 1/float(re.search("SI.hRoiManager.linePeriod = (\d+.?\d*)", md).groups()[0])
    props['pixel_rate'] = 1/float(re.search("SI.hScan2D.scanPixelTimeMean = (\d+.?\d*e-\d+)", md).groups()[0])
    props['bidi_phase_offset'] = float(re.search("SI.hScan2D.linePhase = (\d+.?\d*e-\d+)", md).groups()[0])
    props['bin_factor'] = int(re.search("SI.hScan2D.pixelBinFactor = (\d+)", md).groups()[0])
    
    props['n_channels'] = len(re.search("SI.hChannels.channelSave = \[((?:\d*;?)+)\]", md).groups()[0].split(';'))
    props['frame_shape'] =(int(re.search("SI.hRoiManager.linesPerFrame = (\d+)", md).groups()[0]), int(re.search("SI.hRoiManager.pixelsPerLine = (\d+)", md).groups()[0]))

    #TODO: check this for non-square imaging fields
    um_scanfield = re.search("SI.hRoiManager.imagingFovUm = \[(-?\d+.?\d*) (-?\d+.?\d*);-?\d+.?\d* -?\d+.?\d*;(-?\d+.?\d*) (-?\d+.?\d*);-?\d+.?\d* -?\d+.?\d*\]", md).groups()
    props['microns_per_pixel'] = abs(float(um_scanfield[3]) - float(um_scanfield[1])) / props['frame_shape'][0],abs(float(um_scanfield[2]) - float(um_scanfield[0])) / props['frame_shape'][1]
    
    t = [float(x) for x in re.search("epoch = \[((?:\d+.?\d+,?)+)\]", tif.description(0)).groups()[0].split(',')]
    t.append(t[-1] % 1 * 1e6)
    t = [int(ti) for ti in t]
    offset = 10 * (datetime.datetime(*t) - DOTNET_EPOCH) // US

    props['time_stamps'] = np.empty(len(tif)//4, dtype=np.uint64)
    for i in np.arange(0, len(tif), props['n_channels']):
        props['time_stamps'][i//props['n_channels']] = int(float(re.search("frameTimestamps_sec = (\d+.?\d+)", tif.description(i)).groups()[0]) * 1e7) + offset
    
    data = tif.data().reshape(-1,props['n_channels'],*props['frame_shape']).astype(np.int16)

    if raw_binary_path is not None:
        if not os.path.exists(os.path.join(raw_binary_path, 'functional.bin')):
            func = suite2p.io.BinaryRWFile(Lx=props['frame_shape'][1], Ly=props['frame_shape'][0], filename=os.path.join(raw_binary_path, 'functional.bin'))
            func.write(data[:,functional_channel,:,:])
        if not os.path.exists(os.path.join(raw_binary_path, 'anatomy.bin')):
            anat = suite2p.io.BinaryRWFile(Lx=props['frame_shape'][1], Ly=props['frame_shape'][0], filename=os.path.join(raw_binary_path, 'anatomy.bin'))
            anat.write(data[:,anatomy_channel,:,:])
    
    return data, props

def get_trigger_times(data, props, tracker_channel: int=3, projector_frame_rate: float=60.0):

    trigger = data[:,tracker_channel,:,:].copy() #TODO: should this be a copy?
    trigger[:,1::2,:] = np.flip(trigger[:,1::2,:], axis=2) #bidi scanning
    trigger = trigger.flatten()
    ## Learn the trigger levels and noise using k-means clustering on the trigger channel
    
    w = 4 #window width
    a = trigger[w:] - trigger[:-w]

    IFI = props['pixel_rate'] /  props['frame_rate'] - props['frame_shape'][0] * props['frame_shape'][1]
    ILI = props['pixel_rate'] /  props['line_rate'] -  props['frame_shape'][1]
    mindist = (
        props['pixel_rate'] / projector_frame_rate
        - ILI * np.ceil(props['line_rate'] / projector_frame_rate)
        - IFI * np.ceil(props['frame_rate'] / projector_frame_rate)
        )

    flips,_ = find_peaks(np.abs(a), distance=mindist*.95, prominence=300, wlen=16*600 / props['bin_factor'])
    flips = flips - w//2
    
    flips = np.insert(flips,-1,len(trigger)-1)

    return trigger,flips#.astype(np.uint64)

def cluster_trigger_events(trigger, flips, projector_frame_rate: float=60.0):   

    l = np.empty((len(flips),1))

    l[0] = trigger[:flips[0]][-int(min(flips[0], props['pixel_rate']//projector_frame_rate)):].mean()
    for i in range(1,len(flips)-1):
        l[i] = trigger[(flips[i-1]):flips[i]].mean() # average level of the signal between flips
    l[-1] = trigger[flips[-1]:][:int(min(flips[-1]-flips[-2], props['pixel_rate']//projector_frame_rate))].mean()
    

    km = KMeans(n_clusters=4, algorithm='elkan').fit(l) #we assume we can't learn the start trigger
    fl = km.predict(l)

    # Reorder the predictions in ascending voltage order
    idx = np.argsort(km.cluster_centers_.flatten())
    lut = np.zeros_like(idx)
    lut[idx] = np.arange(len(idx))
    fl = lut[fl]


    ## Calculate the frame number from the trigger labels
    # flips = np.where(np.diff(fl))[0]

    # fl = np.insert(fl,(-1),(0))    

    #flips : the indices of frame boundaries
    #fl    : the type of frame, [0,4)
    #fi    : for each sample, the frame count since the epoch onset

    fi = np.empty_like(trigger, dtype=int)
    fi[:flips[0]] = -1

    fli = np.empty_like(trigger, dtype=int)
    fli[:flips[0]] = -1


    rm = []
    last_flip = flips[0]
    last_fl = 0
    for i in range(1,len(fl)):
        if (flips[i] - flips[i-1] < props['pixel_rate']/(projector_frame_rate * 7/3)):
            #false transition: frames cannot come this close together
            rm.append(i) # keep track of this frame to delete later
            continue
        # elif fl[flips[i]] == 4: # this is the start of a symphony epoch
        #     fi[last_flip:flips[i]] = 3 # cannot distinguish the first 4 frames

        elif fl[i] == 0: # frame is divisible by 4
            fi[last_flip:flips[i]] = ((fi[last_flip-1] + 4)//4)*4
            fli[last_flip:flips[i]] = last_flip
        elif fl[i] == 3: # frame % 4 == 1
            if fl[last_fl] == 0 and flips[i] - flips[last_fl-1] > props['pixel_rate']/projector_frame_rate * 5: #either a skip, or the start trigger
                fi[last_flip:flips[i]] = 3 # cannot distinguish the first 4 frames
                fli[last_flip:flips[i]] = last_flip
            else:
                fi[last_flip:flips[i]] = ((fi[last_flip-1] + 3)//4)*4 + 1
                fli[last_flip:flips[i]] = last_flip
        elif fl[i] == 1: # frame % 4 == 2
            fi[last_flip:flips[i]] = ((fi[last_flip-1] + 2)//4)*4 + 2
            fli[last_flip:flips[i]] = last_flip
        elif fl[i] == 2: # frame % 4 == 3
            fi[last_flip:flips[i]] = ((fi[last_flip-1] + 1)//4)*4 + 3
            fli[last_flip:flips[i]] = last_flip
        
        last_flip = flips[i]
        last_fl = i

    fi[-1] = fi[-2]
    fli[-1] = fli[-2]

    props['flips'] = np.delete(flips, rm)
    props['frame_type'] = np.delete(fl, rm)
    props['frame_index'] = fi
    props['flip_index'] = fli

    props['dropped'] = np.where(np.diff(fi)>1)[0]

    props['dropped'] = props['flips'][np.diff(fi[props['flips']], append=0) > 1]
    
    #TODO: known issues with dropping:
    #   - doesn't correctly handle consecutive frame drops ?
    #   - we frequently drop the second-to-last frame?
    #   - the onset of the first trial increments from -1 to 1 ?
    props['dropped'] = props['dropped'][fi[props['dropped']-1] >= 0]
    #   - the offset of the first frame of each trial increments by 3
    props['dropped'] = props['dropped'][~np.logical_and(fi[props['dropped']] == 0,fi[props['dropped']+1] == 3)]

    props['trials'] = props['flips'][np.diff(fi[props['flips']], prepend=-1) < 0]
    if fi[props['flips'][0]-1] == -1 and fi[props['flips'][0]+1] == 0:
        np.insert(props['trials'], 0, props['flips'][0])


    if len(props['dropped']):
        print(f"Dropped {len(props['dropped'])} frames in {len(props['trials'])} trials ({len(props['dropped'])/(len(props['flips']) + len(props['dropped'])) * 100:.03f}%)")
    return props

def align_frames_to_epochs(epochs: pd.DataFrame, props: Dict[str, Any]) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    '''Align trials between symphony and scanimage
    '''
    
    epochs = epochs.copy()
    props = props.copy()

    time_delta = np.abs(np.expand_dims(epochs['start_time'],axis=0) - np.expand_dims(props['time_stamps'][props['trials']//(props['frame_shape'][0] * props['frame_shape'][1])],axis=1))
    idx_e = time_delta.argmin(axis=0)
    idx_i = time_delta.argmin(axis=1)

    #TODO: isnt this supposed to be argwhere??
    idx_i_missing = np.where(np.take_along_axis(time_delta, np.expand_dims(idx_i, axis=1), axis=1) > 1e7)[1]
    idx_i = np.delete(idx_i, idx_i_missing)

    idx_e_missing = np.where(np.take_along_axis(time_delta, np.expand_dims(idx_e, axis=0), axis=0) > 1e7)[1]
    idx_e = np.delete(idx_e, idx_e_missing)
    idx_e_idx = idx_e.argsort()

    epochs = epochs.drop(index = idx_e_missing).reset_index(drop=True).reindex(idx_e_idx).reset_index(drop=True)

    props['trials']= np.delete(props['trials'], idx_i_missing)

    return epochs, props
    
def time_events(epochs: pd.DataFrame, props: Dict[str, Any]) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, ]:
    #TODO: don't like this... should be a class along with symphony parser
    # and should return dataframe

    if len(np.unique(epochs['spotPreFrames'])) > 1 \
        or len(np.unique(epochs['spotStimFrames'])) > 1 \
        or len(np.unique(epochs['spotTailFrames'])) > 1:
        raise Exception('Spots of different durations are present in the recording')

    spot_trial = int(epochs[['spotPreFrames','spotStimFrames','spotTailFrames']].iloc[0].sum())

    spots = {
        'cx': np.array([]),
        'cy': np.array([]),
        'stim_on': np.array([], dtype=int),
        'stim_on_i': np.array([], dtype=int),
        'stim_off': np.array([], dtype=int),
        'stim_off_i': np.array([], dtype=int),
        'pre_on': np.array([], dtype=int),
        'pre_on_i': np.array([], dtype=int),
        'tail_off': np.array([], dtype=int),
        'tail_off_i': np.array([], dtype=int),
        'intensity': np.array([]),
        'led_level': np.array([], dtype=int),
        'ndf': np.array([], dtype=int)
        }
        
    chirps = {
        'stim_on': np.array([], dtype=int),
        'stim_on_i': np.array([], dtype=int),
        'stim_off': np.array([], dtype=int),
        'stim_off_i': np.array([], dtype=int),
        'pre_on': np.array([], dtype=int),
        'pre_on_i': np.array([], dtype=int),
        'tail_off': np.array([], dtype=int),
        'tail_off_i': np.array([], dtype=int),
        'intensity': np.array([]),
        'led_level': np.array([], dtype=int),
        'ndf': np.array([], dtype=int)
        }
        
    bars = {
        'theta': np.array([]),
        'stim_on': np.array([], dtype=int),
        'stim_on_i': np.array([], dtype=int),
        'stim_off': np.array([], dtype=int),
        'stim_off_i': np.array([], dtype=int),
        'pre_on': np.array([], dtype=int),
        'pre_on_i': np.array([], dtype=int),
        'tail_off': np.array([], dtype=int),
        'tail_off_i': np.array([], dtype=int),
        'intensity': np.array([]),
        'led_level': np.array([], dtype=int),
        'ndf': np.array([], dtype=int)
        }


    flips = props['flips']
    trials = props['trials']
    fi = props['frame_index']

    # for i, (e, tr) in enumerate(zip(epochs, trials)):
    for (i,e), tr in zip(epochs.iterrows(),trials):
        #get the index of the trial start in the frame toggle array
        flip_i =np.where(flips==tr)[0][0]
        if e['trial_type']==b'field':
            pre_i = 0
            on_i = int(epochs['spotPreFrames'][0])
            off_i = int(epochs['spotStimFrames'][0]  + on_i)
            # tail_i = 90
            j = 0
            while (i < len(trials)-1 and flips[flip_i + off_i] < trials[i+1]) or (i == len(trials)-1 and flip_i + off_i < len(flips)):
                #the current index will be the beginning of the 'preTime'
                spots['pre_on'] = np.append(spots['pre_on'], flips[flip_i + np.argwhere(fi[flips[flip_i:]]>=pre_i)[0,0]])
                spots['pre_on_i'] = np.append(spots['pre_on_i'], fi[spots['pre_on'][-1]])
                
                #find the index where the spot goes on (15 frames) and off (30 frames)
                spots['stim_on'] = np.append(spots['stim_on'], flips[flip_i + np.argwhere(fi[flips[flip_i:]]>=on_i)[0,0]])
                spots['stim_on_i'] = np.append(spots['stim_on_i'], fi[spots['stim_on'][-1]])
                
                spots['stim_off'] = np.append(spots['stim_off'], flips[flip_i + np.argwhere(fi[flips[flip_i:]]>=off_i)[0,0]])
                spots['stim_off_i'] = np.append(spots['stim_off_i'], fi[spots['stim_off'][-1]])
                

                spots['cx'] = np.append(spots['cx'], e['cx'][j])
                spots['cy'] = np.append(spots['cy'], e['cy'][j])

                spots['intensity'] = np.append(spots['intensity'], e['spotIntensity'])
                spots['led_level'] = np.append(spots['led_level'], int(e['uvLED']))
                spots['ndf'] = np.append(spots['ndf'], int(e['NDF']))               
            
                #the next trial starts 90 frames later...
                j+=1
                pre_i += spot_trial
                on_i += spot_trial
                off_i += spot_trial

                to = np.argwhere(fi[flips[flip_i:]]>=pre_i)
                if len(to):
                    spots['tail_off'] = np.append(spots['tail_off'], flips[flip_i + to[0,0]])
                else: #the video ends early?
                    spots['tail_off'] = np.append(spots['tail_off'], flips[-1])
                spots['tail_off_i'] = np.append(spots['tail_off_i'], fi[spots['tail_off'][-1]])
                
        elif e['trial_type']==b'chirp':
            on_i = 60*2
            off_i = 60*30
            # tail_i = 90
            
            chirps['pre_on'] = np.append(chirps['pre_on'], flips[flip_i])
            chirps['pre_on_i'] = np.append(chirps['pre_on_i'], fi[chirps['pre_on'][-1]])
                

            chirps['stim_on'] = np.append(chirps['stim_on'], flips[flip_i + np.argwhere(fi[flips[flip_i:]]>=(60*2))[0,0]])
            chirps['stim_on_i'] = np.append(chirps['stim_on_i'], fi[chirps['stim_on'][-1]])

            chirps['stim_off'] = np.append(chirps['stim_off'], flips[flip_i + np.argwhere(fi[flips[flip_i:]]>=(60*30))[0,0]])
            chirps['stim_off_i'] = np.append(chirps['stim_off_i'], fi[chirps['stim_off'][-1]])
        
            tail_off = np.argwhere(fi[flips[flip_i:]]>=(60*35))
            if len(tail_off):
                tail_off = flip_i + tail_off[0,0]
            else:
                tail_off = -1
            chirps['tail_off'] = np.append(chirps['tail_off'], flips[tail_off])
            chirps['tail_off_i'] = np.append(chirps['tail_off_i'], fi[chirps['tail_off'][-1]])

            

            chirps['intensity'] = np.append(chirps['intensity'], e['chirpIntensity'])
            chirps['led_level'] = np.append(chirps['led_level'], int(e['uvLED']))
            chirps['ndf'] = np.append(chirps['ndf'], int(e['NDF']))

        elif e['trial_type']==b'bars':
            pre_i = 0
            on_i = 15
            off_i = 195
            #210

            j = 0
            while (i < len(trials)-1 and flips[flip_i + off_i] < trials[i+1]) or (i == len(trials)-1 and flip_i + off_i < len(flips)):
                #the current index will be the beginning of the 'preTime'
                bars['pre_on'] = np.append(bars['pre_on'], flips[flip_i + np.argwhere(fi[flips[flip_i:]]>=pre_i)[0,0]])
                bars['pre_on_i'] = np.append(bars['pre_on_i'], fi[bars['pre_on'][-1]])
                
                #find the index where the spot goes on (15 frames) and off (30 frames)
                bars['stim_on'] = np.append(bars['stim_on'], flips[flip_i + np.argwhere(fi[flips[flip_i:]]>=on_i)[0,0]])
                bars['stim_on_i'] = np.append(bars['stim_on_i'], fi[bars['stim_on'][-1]])
                
                bars['stim_off'] = np.append(bars['stim_off'], flips[flip_i + np.argwhere(fi[flips[flip_i:]]>=off_i)[0,0]])
                bars['stim_off_i'] = np.append(bars['stim_off_i'], fi[bars['stim_off'][-1]])

                bars['theta'] = np.append(bars['theta'], e['theta'][j])

                bars['intensity'] = np.append(bars['intensity'], e['barIntensity'])
                bars['led_level'] = np.append(bars['led_level'], int(e['uvLED']))
                bars['ndf'] = np.append(bars['ndf'], int(e['NDF']))                
            
                #the next trial starts 90 frames later...
                j+=1
                pre_i += 210
                on_i += 210
                off_i += 210

                to = np.argwhere(fi[flips[flip_i:]]>=pre_i)
                if len(to):
                    bars['tail_off'] = np.append(bars['tail_off'], flips[flip_i + to[0,0]])
                else: #the video ends early?
                    bars['tail_off'] = np.append(bars['tail_off'], flips[-1])
                bars['tail_off_i'] = np.append(bars['tail_off_i'], fi[bars['tail_off'][-1]])

    spot_durations = (np.asarray(spots['stim_off_i'], dtype=int) - np.asarray(spots['stim_on_i'], dtype=int)) #/props['pixel_rate'] * 60
    chirp_durations = (np.asarray(chirps['stim_off_i'], dtype=int) - np.asarray(chirps['stim_on_i'], dtype=int)) #/props['pixel_rate'] * 60
    bar_durations = (np.asarray(bars['stim_off_i'], dtype=int) - np.asarray(bars['stim_on_i'], dtype=int)) #/props['pixel_rate'] * 60

    print(f"Detected {len(spots['pre_on'])} spot epochs with spot durations between {min(spot_durations)} and {max(spot_durations)} frames.")
    print(f"Detected {len(chirps['pre_on'])} chirp epochs with chirp durations between {min(chirp_durations)} and {max(chirp_durations)} frames.")
    print(f"Detected {len(bars['pre_on'])} bar epochs with bar durations between {min(bar_durations)} and {max(bar_durations)} frames.")
    
    return (pd.DataFrame.from_dict(e, orient='columns') for e in (spots,chirps,bars))
    
def do_registration(raw_data_file, stack_file, out_path, props):
    #todo...
    register.register(raw_data_file, stack_file, out_path, props=props, make_gif=True)

def segment(bin_path, props, events, mask_file, dark_file, func_channel = 1, data = None):
    '''
    '''
    # if stim_type is None:
    #     stim_type = np.empty_like(stim_timing.T) * np.nan

    #TODO: this assumes the same data shape as the movie... easy enough to fix
    dark_level = suite2p.io.tiff.ScanImageTiffReader(dark_file).data().reshape(-1,props['n_channels'],*props['frame_shape'])[:,func_channel,:,:].flatten().mean()
    
    if data is None:
        data = suite2p.io.BinaryRWFile(*props['frame_shape'], os.path.join(bin_path,'functional.bin')).data

    
    if type(mask_file) == str:
        mask = suite2p.io.tiff.open_tiff(mask_file, True)[0].pages[0].asarray()
    else:
        mask = mask_file

    reg_x, reg_y = register.get_shift(bin_path, mask.shape, props)
    reg_x = reg_x.flatten()
    reg_y = reg_y.flatten()
    # reg_x and reg_y are in the stack coordinates


    F,Y,X = data.shape

    tt = np.zeros((F,Y,X))

    tt[:,0::2,:] += np.arange(0, X)[None, None, :] / props['pixel_rate']
    tt[:,1::2,:] += np.flip(np.arange(0, X)[None, None, :] / props['pixel_rate'], axis=2) #for bidi scanning
    tt += np.arange(0, Y)[None, :, None] / props['line_rate']
    tt += ((props['time_stamps'] - props['time_stamps'][0]) / 1e7)[:, None, None]

    xt = np.tile(np.arange(0,X)[None,None,:], (F,Y,1))
    yt = np.tile(np.arange(0,Y)[None,:,None], (F,1,X))
    ft = np.tile(np.arange(0,F)[:,None,None], (1,Y,X))

    ##flatten the parameters
    xt = xt.flatten()
    yt = yt.flatten()
    ft = ft.flatten()
    tt = tt.flatten()
    dt = data.flatten()


    ##create a dataframe holding all the observations
    n_obs = (events['tail_off'] - events['pre_on']).sum()
    df_raw = pd.DataFrame().assign(        
        fluor = np.empty(n_obs),                    #the raw fluorescence value
        t = np.empty(n_obs),                        #time since the start of the recording
        trial_t = np.empty(n_obs),                  #time since the start of the trial
        # frame_t = np.empty(n_obs),                  #time since the start of the projector frame
        x = np.empty(n_obs),                        #x position in reference stack coordinates
        y = np.empty(n_obs),                        #y position in reference stack coordinates
        col = np.empty(n_obs, dtype=np.uint16),     #x position in movie coordinates 
        line = np.empty(n_obs, dtype=np.uint16),    #y position in movie coordinates
        frame = np.empty(n_obs, dtype=np.uint16),   #t position in movie coordinates
        trial = np.empty(n_obs, dtype=np.uint16),   #trial number from the `events` dataframe
        roi = np.empty(n_obs, dtype=np.uint16),     #roi number from the mask file
        baseline_fluor = np.empty(n_obs),           #average fluorescence in the period preceding the trial start
        dFoF = np.empty(n_obs),                     #normalized fluorescence
        )

    N = 0
    for epoch in events.itertuples():
        trial_pts = min(epoch.tail_off,len(dt)-1) - epoch.pre_on
        trial_out = np.arange(N, N+trial_pts)
        trial_in = np.arange(epoch.pre_on, epoch.pre_on+trial_pts)

        df_raw.loc[trial_out,'t'] = tt[trial_in]
        df_raw.loc[trial_out,'trial_t'] = tt[trial_in] - tt[epoch.stim_on]
        # df_raw.loc[trial_out, 'frame_t'] = tt[trial_in] - tt[props['flips'][trial_in]] #this isn't even right?

        df_raw.loc[trial_out,'col'] = xt[trial_in]
        df_raw.loc[trial_out,'line'] = yt[trial_in]
        df_raw.loc[trial_out,'frame'] = ft[trial_in]

        df_raw.loc[trial_out,'fluor'] = dt[trial_in] - dark_level

        df_raw.loc[trial_out,'trial'] = epoch.Index     

        # df_raw.loc[trial_out,'x'] = reg_x[trial_in]
        # df_raw.loc[trial_out,'y'] = reg_y[trial_in]
        df_raw.loc[trial_out,'x'] = reg_x[np.ravel_multi_index((ft[trial_in], yt[trial_in], xt[trial_in]), (len(reg_x), *props['frame_shape']))]
        df_raw.loc[trial_out,'y'] = reg_y[np.ravel_multi_index((ft[trial_in], yt[trial_in], xt[trial_in]), (len(reg_y), *props['frame_shape']))]    

        N += trial_pts

    df_raw.set_index('t', inplace=True) #NOTE: this assumes that there are no overlapping pre- and tail- times?

    df_raw.drop(df_raw.index[df_raw.x <= -0.5], inplace=True)
    df_raw.drop(df_raw.index[df_raw.y <= -0.5], inplace=True)
    df_raw.drop(df_raw.index[df_raw.x >= mask.shape[1] - 0.5], inplace=True)
    df_raw.drop(df_raw.index[df_raw.y >= mask.shape[0] - 0.5], inplace=True)

    df_raw.loc[:,'roi'] = mask.flat[np.ravel_multi_index((np.round(df_raw['y']).astype(int),np.round(df_raw['x']).astype(int)), mask.shape)]
    df_raw.drop(df_raw.index[df_raw.roi == 0], inplace=True)

    # df_raw.loc[:,'baseline_fluor'] = df_raw.groupby(['trial','roi']).transform(lambda x: x[x['trial_t']<0]['fluor'].mean())
    
    tmp = df_raw.groupby(['trial','roi']).apply(lambda x: pd.Series({'baseline_fluor':x[x['trial_t']<0]['fluor'].mean()}))
    df_raw.loc[:,'baseline_fluor'] = df_raw[['trial','roi']].join(tmp, on=['trial','roi'])['baseline_fluor']
    
    df_raw.loc[:,'dFoF'] = df_raw.apply(lambda x: (x['fluor'] - x['baseline_fluor'])/x['baseline_fluor'], axis=1)

    return df_raw


if __name__ == "__main__":

    symphony_file = r"C:\Users\zfj\data\120822B\120822B.h5"
    raw_file = r"C:\Users\zfj\data\120822B\region2_00001.tif"
    mask_file = r"C:\Users\zfj\data\training\120822B_region2_34_masks.tif"
    dark_file = r"C:\Users\zfj\data\120822B\region2_shuttered_00001.tif"
    stack_file = r"C:\Users\zfj\data\120822B\region2_stack_00001.tif"
    bin_path =r"C:\Users\zfj\data\120822B\func"
    

    blocks = extract_parameters(symphony_file)
    epochs = flatten_epochs(blocks)

    props = process_raw_file(raw_file, bin_path)
    epochs, props = align_frames_to_epochs(epochs, props)

    spots, chirps, bars = time_events(epochs, props)
    bars_df = segment(bin_path, props, bars, mask_file, dark_file, func_channel = 1)

