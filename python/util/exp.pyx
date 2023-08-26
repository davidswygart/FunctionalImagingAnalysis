import numpy as np
from sklearn.cluster import KMeans
import pandas as pd

from .utils import intervals, chirp_pattern

import cython
cimport numpy as cnp
cimport cython

cnp.import_array()


ctypedef cnp.uint8_t uint8_t
ctypedef cnp.uint64_t uint64_t
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def align_frames_to_epochs(cnp.ndarray[uint64_t, ndim=1] frames, cnp.ndarray[uint64_t, ndim=1] starts, cnp.ndarray[uint64_t, ndim=1] ends):
    '''
    Assumes that all three vectors are in ascending order
    '''
    cdef cnp.ndarray[cnp.int64_t, ndim=2] frame_i = np.ones((starts.shape[0],2), dtype=np.int64) * -1

    cdef int i = 0
    cdef int j = 0

    cdef int ne = len(starts)
    cdef int nf = len(frames)

    while frames[i] > starts[j]: # the early frames might occur after the early epochs
        j += 1

    while j < ne:
        while frames[i] < starts[j]:
            i += 1
            if i == nf:
                return frame_i
        frame_i[j,0] = i - 1

        while frames[i] < ends[j]:
            i += 1
            if i == nf:
                return frame_i
        frame_i[j,1] = i - 1
        j += 1
    return frame_i


def cluster_trigger_events(trigger, flips, begin, end, pixel_rate, projector_frame_rate=60.0):
    # trigger = e.trigger[begin:end]
    flps = flips[np.logical_and(flips >= begin, flips < end)] - begin

    lv = np.asarray([trigger[begin:end][flps[i-1]:flps[i]].mean() for i in range(1, len(flps))]).reshape((-1,1))
    
    #ignore the first 10 flips as (possibly) belonging to the start trigger
    km = KMeans(n_clusters=4, algorithm='elkan', n_init=10).fit(lv[10:])
    fl = km.predict(lv)
    idx = np.argsort(km.cluster_centers_.flatten())
    lut = np.zeros_like(idx)
    lut[idx] = np.arange(len(idx))
    fl = [0, *lut[fl]]

    fi = np.empty(len(flps), dtype=int)
    fi[0] = -1

    rm = []
    last_fl = 0
        
    for i in range(1,len(fl)):
        if (flps[i] - flps[i-1] < pixel_rate/(projector_frame_rate * 7/3)) or (i<len(fl)-1 and fl[i] == fl[i+1]):
            #false transition: frames cannot come this close together
            rm.append(i) # keep track of this frame to delete later
            continue
        elif fl[i] == 0: # frame is divisible by 4
            fi[i] = ((fi[last_fl] + 4)//4)*4
        elif fl[i] == 3: # frame % 4 == 3
            if not last_fl or (fi[last_fl] == -1 and flps[i] - flps[last_fl] > pixel_rate/projector_frame_rate * 5): #either a skip, or the start trigger
                fi[i] = 3 # cannot distinguish the first 4 frames
            else:
                fi[i] = ((fi[last_fl] + 3)//4)*4 + 1
        elif fl[i] == 1: # frame % 4 == 1
            fi[i] = ((fi[last_fl] + 2)//4)*4 + 2
        elif fl[i] == 2: # frame % 4 == 2
            fi[i] = ((fi[last_fl] + 1)//4)*4 + 3
        else:
            print('!')
        
        # last_flip = flips[i]
        last_fl = i    

    flps = np.delete(flps, rm)
    fl = np.delete(fl, rm) # not used?
    fi = np.delete(fi, rm)[1:]

    # dropped = np.where(np.diff(fi)>1)[0]

    trial_i = np.logical_and(fi[:-1]<3, fi[1:]==3)

    dropped_i = np.logical_and(np.diff(fi, prepend=-1) > 1, fi>=0)
    dropped_i = np.logical_and(dropped_i[:-1], ~trial_i)
    dropped = flps[np.argwhere(dropped_i)]

    # #TODO: known issues with dropping:
    # #   - doesn't correctly handle consecutive frame drops ?
    # #   - we frequnetly drop the second-to-last frame ?
    # #   - claims that frames are dropped when there are false triggers


    # trials = flips[np.diff(fi[flips], prepend=-1) < 0]
    trials = flps[np.argwhere(trial_i)]

    if len(dropped):
        print(f"Dropped {len(dropped)} frames in {len(trials)} trials ({len(dropped)/(len(flips) + len(dropped)) * 100:.03f}%)")

    return fi, flps + begin, fl, dropped + begin, trials + begin, km


def time_events(epochs):
    cols = {
        'pre_on': int, 'stim_on': int, 'stim_off': int, 'tail_off': int,
        'pre_on_i': int, 'stim_on_i': int, 'stim_off_i': int, 'tail_off_i': int, 
        'epoch_id': int
        }
    scols = {**cols, 'cx':np.float64, 'cy':np.float64}
    spots = pd.DataFrame(columns=scols.keys()).astype(scols, copy=False, errors='ignore')
    bcols = {**cols, 'theta':np.float64}
    bars = pd.DataFrame(columns=bcols.keys()).astype(bcols, copy=False, errors='ignore')
    chirps = pd.DataFrame(columns=cols.keys()).astype(cols, copy=False, errors='ignore')

    for i,epoch in epochs.iterrows():
        flips = epoch.flips
        if flips is None:
            continue
        fi = epoch.frame_indices

        # flip_i =np.where(flips==tr)[0][0] # hmmm...
        # flip_i = flips[0]
        # flip_i = 0

        if epoch.trial_type==b'field':
            spot_trial = int(epoch[['spotPreFrames','spotStimFrames','spotTailFrames']].sum())
            spots_per_epoch = 35*60//spot_trial

            spots_t = pd.DataFrame(np.zeros((spots_per_epoch,len(scols))), columns=scols).astype(scols, copy=False, errors='ignore')
            spots_t.epoch_id = i
            pre_i = 0
            on_i = int(epoch.spotPreFrames)
            off_i = int(epoch.spotStimFrames + on_i)
            tail_i = spot_trial
            j = 0
            k = 0

            while j < spots_per_epoch:#k < len(flips):#, flip in enumerate(flips):
                # associate the timing events with the first frame flip that exceeds the requisite frame count
                # specifically, the timing of the named events for each trial will be linked with a linear index into the movie

                #pre_on, pre_on_i
                pre_on = np.argwhere(fi[k:]>=pre_i)[0,0]
                spots_t.iat[j, 0] = flips[k + pre_on]
                spots_t.iat[j, 4] = fi[k + pre_on]
                
                #stim_on, stim_on_i
                stim_on = np.argwhere(fi[k:]>=on_i)[0,0]
                spots_t.iat[j, 1] = flips[k + stim_on]
                spots_t.iat[j, 5] = fi[k + stim_on]

                #stim_off, stim_off_i
                stim_off = np.argwhere(fi[k:]>=off_i)[0,0]
                spots_t.iat[j, 2] = flips[k + stim_off]
                spots_t.iat[j, 6] = fi[k + stim_off - 1]

                tail_off = np.argwhere(fi[k:]>=tail_i)
                if len(tail_off):
                    tail_off = tail_off[0,0]
                else: #we're missing the last flip
                    tail_off = len(flips) - k - 1
                
                spots_t.iat[j, 3] = flips[k + tail_off]
                spots_t.iat[j, 7] = fi[k + tail_off - 1]
                

                spots_t.iat[j, 9] = epoch.cx[j]
                spots_t.iat[j, 10] = epoch.cy[j]


                #increment the trial counter
                k = tail_i
                pre_i += spot_trial
                on_i += spot_trial
                off_i += spot_trial
                tail_i += spot_trial
            
                j+=1
            
            spots = pd.concat([spots, spots_t])
                
        elif epoch.trial_type==b'chirp':
            chirps_t = pd.DataFrame(np.zeros((1,len(cols))), columns=cols).astype(cols, copy=False, errors='ignore')
            chirps_t.epoch_id = i
            
            chirps_t.pre_on= flips[0]
            chirps_t.pre_on_i= fi[0]

            stim_on = np.argwhere(fi>=60*2)[0,0] #TODO: should be frame rate
            chirps_t.stim_on = flips[stim_on]
            chirps_t.stim_on_i = fi[stim_on]

            stim_off = np.argwhere(fi>=60*30)[0,0]
            chirps_t.stim_off = flips[stim_off]
            chirps_t.stim_off_i = fi[stim_off - 1]
            
            tail_off = np.argwhere(fi>=60*35)
            if len(tail_off):
                tail_off = tail_off[0,0]
            else:
                tail_off = len(flips) - 1
            
            chirps_t.tail_off = flips[tail_off]
            chirps_t.tail_off_i = fi[tail_off - 1]

            
            chirps = pd.concat([chirps, chirps_t])

        elif epoch.trial_type==b'bars':            
            bars_t = pd.DataFrame(np.zeros((10,len(bcols))), columns=bcols).astype(bcols, copy=False, errors='ignore')
            bars_t.epoch_id = i            
            
            pre_i = 0
            on_i = 15
            off_i = 195
            tail_i = 210
            j = 0
            k = 0

            while j < 10:#k < len(flips):#, flip in enumerate(flips):
                # associate the timing events with the first frame flip that exceeds the requisite frame count
                # specifically, the timing of the named events for each trial will be linked with a linear index into the movie

                #pre_on, pre_on_i
                pre_on = np.argwhere(fi[k:]>=pre_i)[0,0]
                bars_t.iat[j, 0] = flips[k + pre_on]
                bars_t.iat[j, 4] = fi[k + pre_on]
                
                #stim_on, stim_on_i
                stim_on = np.argwhere(fi[k:]>=on_i)[0,0]
                bars_t.iat[j, 1] = flips[k + stim_on]
                bars_t.iat[j, 5] = fi[k + stim_on]

                #stim_off, stim_off_i
                stim_off = np.argwhere(fi[k:]>=off_i)[0,0]
                bars_t.iat[j, 2] = flips[k + stim_off]
                bars_t.iat[j, 6] = fi[k + stim_off - 1]

                tail_off = np.argwhere(fi[k:]>=tail_i)
                if len(tail_off):
                    tail_off = tail_off[0,0]
                else: #we're missing the last flip
                    tail_off = len(flips) - k - 1
                
                bars_t.iat[j, 3] = flips[k + tail_off]
                bars_t.iat[j, 7] = fi[k + tail_off - 1]
                

                bars_t.iat[j, 9] = epoch.theta[j]


                #increment the trial counter
                k = tail_i
                pre_i += 210
                on_i += 210
                off_i += 210
                tail_i += 210
            
                j+=1
            
            bars = pd.concat([bars, bars_t])
    return spots, chirps, bars



cdef extern from "mask.cpp":
    void mask(
        double* y,
        double* x,
        double* w,

        int width,
        int n_trials,

        uint64_t* frames,
        uint64_t* n_frames,
        uint64_t* flips,
        uint64_t* trial_type,

        uint64_t* pre_frames,
        uint64_t* stim_frames,
        uint64_t* tail_frames,

        double* frame_rate,
        double* intensity
    )
    # void chirp(uint8_t* pattern, uint8_t intensity, uint64_t frame_rate) {

    # void chirp(
    #     uint8_t* pattern, uint8_t intensity, uint64_t frame_rate
    # )
    void chirp(
        uint8_t* pattern, double intensity, uint64_t frame_rate
    )

def get_chirp(intensity, frame_rate, pattern=None):

    if pattern is None:
        pattern = np.ascontiguousarray(np.empty(2098, dtype=np.uint8))
    else:
        pattern = np.ascontiguousarray(pattern)
    cdef uint8_t[::1] pattern_mv = pattern
    chirp(&pattern_mv[0], intensity, frame_rate)

    return pattern


def mask_artifact(func, events, props, artifact):
    F,Y,X = func.shape
    func[:,1::2,:] = np.flip(func[:,1::2,:], axis=2)
    func = func.flatten().astype(np.double)

    tt = np.zeros((F,Y,X), dtype=np.double)
    tt += np.arange(0, X)[None, None, :] / props['pixel_rate']
    tt += np.arange(0, Y)[None, :, None] / props['line_rate']
    tt += ((props['time_stamps'] - props['time_stamps'][0]) / 1e7)[:, None, None]
    tt = tt.flatten()

    

    # events[['trial_type','intensity', 'spotPreFrames','spotStimFrames','frameRate']].to_numpy()
    #flatten the frame_indices array
    fi = np.concatenate(events.frame_indices.tolist())
    nfi = events.apply(lambda x : len(x.frame_indices), axis=1).to_numpy()
    flps = np.concatenate(events.flips.tolist())
    trial_type = (1*(events['trial_type']==b'field') + 2*(events['trial_type']==b'bars')).to_numpy()

    frame_rate = events['frameRate'].to_numpy().astype(np.double)

    intensity = events['intensity'].to_numpy().copy()
    intensity[np.logical_not(trial_type)] *= 2
    # intensity = (intensity * 255).astype(np.uint8)
    spot_pre = events['spotPreFrames'].to_numpy().astype(np.double)
    spot_stim = events['spotStimFrames'].to_numpy().astype(np.double)
    spot_tail = events['spotTailFrames'].to_numpy().astype(np.double)

    

    return _mask_artifact(func, 
                   tt, 
                   artifact, 
                   props['frame_shape'][1], 
                   fi.astype(np.uint64), 
                   nfi.astype(np.uint64), 
                   flps.astype(np.uint64), 
                   trial_type.astype(np.uint64), 
                   frame_rate.astype(np.double), 
                   intensity.astype(np.double), 
                   spot_pre.astype(np.uint64), 
                   spot_stim.astype(np.uint64),
                   spot_tail.astype(np.uint64),
                   )

def _mask_artifact(
        func,
        tt,
        artifact,
        image_width,
        frame_indices,
        nframes,
        flip_indices,
        trial_type,
        frame_rate,
        intensity,
        spotPreFrames,
        spotStimFrames,
        spotTailFrames,
):
    
    
    # if not func.flags['C_CONTINUOUS']:
    func = np.ascontiguousarray(func)
    cdef double[::1] func_mv = func

    # if not tt.flags['C_CONTINUOUS']:
    tt = np.ascontiguousarray(tt)
    cdef double[::1] tt_mv = tt      
    
    # if not artifact.flags['C_CONTINUOUS']:
    artifact = np.ascontiguousarray(artifact)
    cdef double[::1] artifact_mv = artifact.flatten()
    
    
    # if not frame_rate.flags['C_CONTINUOUS']:
    frame_rate = np.ascontiguousarray(frame_rate)
    cdef double[::1] frame_rate_mv = frame_rate


    # if not intensity.flags['C_CONTINUOUS']:
    intensity = np.ascontiguousarray(intensity)
    cdef double[::1] intensity_mv = intensity



    # if not frame_indices.flags['C_CONTINUOUS']:
    frame_indices = np.ascontiguousarray(frame_indices)
    cdef uint64_t[::1] frame_indices_mv = frame_indices

    # if not nframes.flags['C_CONTINUOUS']:
    nframes = np.ascontiguousarray(nframes)
    cdef uint64_t[::1] nframes_mv = nframes

    # if not flip_indices.flags['C_CONTINUOUS']:
    flip_indices = np.ascontiguousarray(flip_indices)
    cdef uint64_t[::1] flip_indices_mv = flip_indices

    # if not trial_type.flags['C_CONTINUOUS']:
    trial_type = np.ascontiguousarray(trial_type)
    cdef uint64_t[::1] trial_type_mv = trial_type

    # if not spotPreFrames.flags['C_CONTINUOUS']:
    spotPreFrames = np.ascontiguousarray(spotPreFrames)
    cdef uint64_t[::1] spotPreFrames_mv = spotPreFrames

    # if not spotStimFrames.flags['C_CONTINUOUS']:
    spotStimFrames = np.ascontiguousarray(spotStimFrames)
    cdef uint64_t[::1] spotStimFrames_mv = spotStimFrames

    
    # if not spotTailFrames.flags['C_CONTINUOUS']:
    spotTailFrames = np.ascontiguousarray(spotTailFrames)
    cdef uint64_t[::1] spotTailFrames_mv = spotTailFrames


    print("Got to c call")

    mask(
        &func_mv[0],
        &tt_mv[0],
        &artifact_mv[0],

        image_width,
        len(nframes),

        &frame_indices_mv[0],
        &nframes_mv[0],
        &flip_indices_mv[0],
        &trial_type_mv[0],

        &spotPreFrames_mv[0],
        &spotStimFrames_mv[0],
        &spotTailFrames_mv[0],

        &frame_rate_mv[0],
        &intensity_mv[0],
    )

    return func




## func = func.copy().astype(np.double)

## ctypedef cnp.uint64_t uint64_t

# ctypedef cnp.float64_t double_t

# @cython.boundscheck(False) # turn off bounds-checking for entire function
# @cython.wraparound(False)  # turn off negative index wrapping for entire function
# def _mask_artifact(
#     cnp.ndarray[double_t, ndim=1] func, 
#     cnp.ndarray[double_t, ndim=1] tt, 
#     cnp.ndarray[double_t, ndim=2] lookup, 
#     int image_width, 
#     cnp.ndarray[uint64_t, ndim=1] frame_indices,
#     cnp.ndarray[uint64_t, ndim=1] nframes,  
#     cnp.ndarray[uint64_t, ndim=1] flip_indices,
#     cnp.ndarray[uint64_t, ndim=1] trial_type,  
#     cnp.ndarray[double_t, ndim=1] frame_rate, 
#     cnp.ndarray[double_t, ndim=1] intensity, 
#     cnp.ndarray[uint64_t, ndim=1] spotPreFrames,  
#     cnp.ndarray[uint64_t, ndim=1] spotPreStimFrames,  
    
#     ):

#     # cdef double_t lasta = np.nan
#     cdef int i, j, k, l, a, b, lasta

#     cdef cnp.ndarray[cnp.uint8_t, ndim=1] pattern = np.empty(2200, dtype=np.uint8)


#     k = 0
#     l = 0

#     from tqdm import tqdm
#     for j in tqdm(range(len(nframes))):
#         # print(j)
#         # row = events.iloc[j]
#         lasta = frame_indices[l]


#         if trial_type[j]==0: # chirp
#             pattern[:nframes[j]] = chirp_pattern(intensity[j] * 2, frame_rate[j])[frame_indices[l:(l+nframes[j])]]

#         elif trial_type[j]==1: # spots
#             pattern[:] = 0
#             pattern[(frame_indices[l:(l+nframes[j])] >= spotPreFrames[j]) & (frame_indices[l:(l+nframes[j])] < (spotPreStimFrames[j]))] = int(255*intensity[j])
        
#         elif trial_type[j]==2: # bars
#             pattern[:] = 0
#             pattern[np.logical_and(frame_indices[l:(l+nframes[j])] >= 15, frame_indices[l:(l+nframes[j])] < 195)] = int(255*intensity[j])
        

#         for i,(a,b) in enumerate(intervals(flip_indices[(l+j):(l+nframes[j]+j+1)])): #flip indices has one more element than frame_indices per epoch...
#             if i>=2098:
#                 break
            
#             if (not (a - 1) % image_width) or (not (b - 1) % image_width):
                
#                 func[a:b] = np.nan
#                 # print('.',end='')
#                 lasta = a
#                 continue

#             lu = lookup[pattern][i]
#             lu[lu==0] = np.nan
            
#             inds = ((tt[a:b] - tt[a]) / (tt[b] - tt[a]) * 256).astype(int)
#             for ind,ab in zip(inds, np.arange(a,b)):
#                 # f = func[ab]
#                 # spat[i,j,ind] += f
#                 # npat[i,j,ind] += 1
#                 func[ab] *= lu[ind]

#             if not np.any(inds == 1) and not np.any(inds == 8):
#             # if npat[i,j,1] == 0 and npat[i,j,8] == 0:
#                 # print('!', end='')
#                 # spat[i,j,:] = np.nan
#                 # npat[i,j,:] = np.nan

#                 # spat[i-1,j,:] = np.nan
#                 # npat[i-1,j,:] = np.nan

#                 func[a:b] = np.nan
#                 func[lasta:a] = np.nan

#             # nmask = (lookup[pattern[:2098]]).astype(float)
#             # nmask[nmask==0] = np.nan
#             # spat2[i,j,:] = spat[i,j,:]
#             # spat2[i,j,:] *= lu
#             lasta = a
#         k+=1
#         l+= nframes[j]