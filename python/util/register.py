import suite2p

import matplotlib.pyplot as plt
from matplotlib import animation

import re
import os
from skimage.transform import resize
import numpy as np

import pickle
from tqdm import tqdm

def compute_reference_masks(refImg, ops=suite2p.default_ops()):
    '''
    Override of built-in, which doesn't properly forward ops for z stacks
    '''
    ### ------------- compute registration masks ----------------- ###
    if isinstance(refImg, list):
        refAndMasks_all = []
        for rimg in refImg:
            refAndMasks = compute_reference_masks(rimg, ops) #<- bug fix from built-in
            refAndMasks_all.append(refAndMasks)
        return refAndMasks_all
    else:
        maskMul, maskOffset = suite2p.registration.rigid.compute_masks(
            refImg=refImg,
            maskSlope=ops['spatial_taper'] if ops['1Preg'] else 3 * ops['smooth_sigma'],
        )
        cfRefImg = suite2p.registration.rigid.phasecorr_reference(
            refImg=refImg,
            smooth_sigma=ops['smooth_sigma'],
        )
        Ly, Lx = refImg.shape
        if ops.get('nonrigid'):
            blocks = suite2p.registration.nonrigid.make_blocks(Ly=Ly, Lx=Lx, block_size=ops['block_size'])

            maskMulNR, maskOffsetNR, cfRefImgNR = suite2p.registration.nonrigid.phasecorr_reference(
                refImg0=refImg,
                maskSlope=ops['spatial_taper'] if ops['1Preg'] else 3 * ops['smooth_sigma'], # slope of taper mask at the edges
                smooth_sigma=ops['smooth_sigma'],
                yblock=blocks[0],
                xblock=blocks[1],
            )
        else:
            maskMulNR, maskOffsetNR, cfRefImgNR, blocks = [], [], [], []

        return maskMul, maskOffset, cfRefImg, maskMulNR, maskOffsetNR, cfRefImgNR, blocks

def registration_wrapper(func_in, func_out, anat_in, anat_out, ref_stack, ops=suite2p.default_ops()):
    stack_shape = ref_stack.shape[1:]
    data_shape = anat_in.shape[1:]
    ds_factor = (stack_shape[0]//data_shape[0], stack_shape[1]//data_shape[1])

    nr,rmin,rmax = suite2p.registration.register.normalize_reference_image([plane for plane in ref_stack])
    
    refAndMasks = compute_reference_masks(nr, ops=ops)
    maskMul, maskOffset, cfRefImg, maskMulNR, maskOffsetNR, cfRefImgNR, blocks = list(map(list, zip(*refAndMasks)))
    blocks_ds = ([y // ds_factor[0] for y in blocks[0][0]], [x // ds_factor[1] for x in blocks[0][1]])
    n_blocks = np.product(blocks[0][2])

    frames = np.empty((ops['batch_size'], *stack_shape)) #this will be discarded each iteration and we will redo...
    
    yoff_r = np.empty(anat_in.shape[0], dtype=np.float32)
    xoff_r = np.empty(anat_in.shape[0], dtype=np.float32)
    cmax_r = np.empty(anat_in.shape[0])
    yoff = np.empty((anat_in.shape[0], n_blocks), dtype=np.float32)
    xoff = np.empty((anat_in.shape[0], n_blocks), dtype=np.float32)
    cmax = np.empty((anat_in.shape[0], n_blocks))
    zoff = np.empty(anat_in.shape[0], dtype=int)
    zcmax = np.empty((anat_in.shape[0], len(ref_stack)))
    

    for i in tqdm(np.arange(0, anat_in.shape[0], ops['batch_size'])):
        n_frames = min(ops['batch_size'], anat_in.shape[0] - i)

        batch = resize(anat_in[i:i+n_frames],(n_frames,*stack_shape), preserve_range=True).astype(np.int16)
        # frames[:n_frames],
        # yoff_r[i:i+n_frames],
        # xoff_r[i:i+n_frames],
        # cmax_r[i:i+n_frames], 
        # yoff[i:i+n_frames,:], 
        # xoff[i:i+n_frames,:], 
        # cmax[i:i+n_frames,:],
        # (zoff[i:i+n_frames], 
        # zcmax[i:i+n_frames,:]) = suite2p.registration.register.register_frames(
        #     outputs, batch, rmin=rmin, rmax=rmax, nZ=len(ref_stack), ops=ops
        # )
        outputs = suite2p.registration.register.register_frames(
            refAndMasks, batch, rmin=rmin, rmax=rmax, nZ=len(ref_stack), ops=ops
        )
        frames[:n_frames],yoff_r[i:i+n_frames],xoff_r[i:i+n_frames],cmax_r[i:i+n_frames],yoff[i:i+n_frames,:],xoff[i:i+n_frames,:],cmax[i:i+n_frames,:],z = outputs
        zoff[i:i+n_frames], zcmax[i:i+n_frames,:] = z

        anat_out.write(suite2p.registration.nonrigid.transform_data( 
            data=anat_in[i:i+n_frames],
            yblock=blocks_ds[0],
            xblock=blocks_ds[1],
            nblocks=blocks[0][2],
            ymax1=(yoff[i:i+n_frames, :] + yoff_r[i:i+n_frames, None]) / ds_factor[0],
            xmax1=(xoff[i:i+n_frames, :] + xoff_r[i:i+n_frames, None]) / ds_factor[1],
        ))

        if func_in is not None and func_out is not None:
            func_out.write(suite2p.registration.nonrigid.transform_data( 
                data=func_in[i:i+n_frames],
                yblock=blocks_ds[0],
                xblock=blocks_ds[1],
                nblocks=blocks[0][2],
                ymax1=(yoff[i:i+n_frames, :] + yoff_r[i:i+n_frames, None]) / ds_factor[0],
                xmax1=(xoff[i:i+n_frames, :] + xoff_r[i:i+n_frames, None]) / ds_factor[1],
            ))
    return yoff_r, xoff_r, cmax_r, yoff, xoff, cmax, zoff, zcmax

def register(movie_path, stack_path, out_path,  props, anatomy_chan=0, func_chan=1, make_gif=True, overwrite_raw=True):

    stack,_ = suite2p.io.tiff.open_tiff(stack_path,sktiff=0)
    
    n_frames_per_plane = int(re.search("SI.hStackManager.framesPerSlice = (\d+)", stack.metadata()).groups()[0]) // int(re.search("SI.hScan2D.logAverageFactor = (\d+)", stack.metadata()).groups()[0])
    stack_size =(int(re.search("SI.hRoiManager.linesPerFrame = (\d+)", stack.metadata()).groups()[0]), int(re.search("SI.hRoiManager.pixelsPerLine = (\d+)", stack.metadata()).groups()[0]))
    stack_n_chans = len(re.search("SI.hChannels.channelSave = \[((?:\d+;?)+)\]", stack.metadata()).groups()[0].split(';'))

    stack_data = stack.data().reshape((-1,n_frames_per_plane,stack_n_chans,*stack_size))[:,:,anatomy_chan].mean(axis=1)

    in_path, in_file = os.path.split(movie_path)

    if os.path.exists(os.path.join(out_path, 'functional.bin')):
        raw_file = os.path.join(out_path, 'functional.bin')
        reg_file = os.path.join(out_path, 'functional_reg.bin')
        do_func = True
    else:
        raw_file = os.path.join(out_path, 'anatomy.bin')
        reg_file = os.path.join(out_path, 'anatomy_reg.bin')
        do_func = False


    ops = suite2p.default_ops()
    ops.update(dict(
        nplanes = 1, #will this work?
        tiff_list = [in_file],
        data_path = [in_path],
        raw_file = raw_file,
        raw_file_chan2 = os.path.join(out_path, 'anatomy.bin'),   
        reg_file = reg_file,
        reg_file_chan2 = os.path.join(out_path, 'anatomy_reg.bin'),   
        save_path0 = out_path,#os.path.join(out_path,'movie.bin'),
        save_folder = out_path,
        nchannels = props['n_channels'],
        align_by_chan = anatomy_chan + 1,
        functional_chan = func_chan + 1,
        nonrigid = True,
        keep_movie_raw = True,
        input_format = 'tif'
    ))

    if do_func:
        func = suite2p.io.BinaryRWFile(Lx=props['frame_shape'][1], Ly=props['frame_shape'][0], filename=ops['raw_file'])        
        reg_f = suite2p.io.BinaryRWFile(Lx=props['frame_shape'][1], Ly=props['frame_shape'][0], filename=ops['reg_file'])
    else:
        func = None
        reg_f = None
    anat = suite2p.io.BinaryRWFile(Lx=props['frame_shape'][1], Ly=props['frame_shape'][0], filename=ops['raw_file_chan2'])
    reg_a = suite2p.io.BinaryRWFile(Lx=props['frame_shape'][1], Ly=props['frame_shape'][0], filename=ops['reg_file_chan2'])
    
    
    outputs = registration_wrapper(func, reg_f, anat, reg_a, stack_data, ops=ops)
    yoff_r, xoff_r, cmax_r, yoff, xoff, cmax, zoff, zcmax = outputs
    reg_data = {
        'y_rigid':yoff_r,
        'x_rigid':xoff_r,
        'y_nr': yoff,
        'x_nr': xoff,
        'cmax_rigid': cmax_r,
        'cmax_nr' : cmax,
        'z':zoff,
        'z_cmax':zcmax,
    }

    # suite2p.registration.register.save_registration_outputs_to_ops(outputs, ops)
    np.save(os.path.join(out_path,'ops.npy'), ops, allow_pickle=True)
    np.save(os.path.join(out_path,'reg.npy'), reg_data, allow_pickle=True)    

    if make_gif:
        fig = plt.figure()
        i = np.random.randint(0,anat.shape[0])
        im = plt.imshow(anat[i].squeeze(), aspect=props['frame_shape'][1]/props['frame_shape'][0]) #assumes nonsquare pix...
        plt.clim(np.percentile(anat[i],10), np.percentile(anat[i],90))


        def animate(i):
            im.set_data(anat[np.random.randint(0,anat.shape[0])].squeeze())
            return [im]

        anim = animation.FuncAnimation(fig, animate, frames = 300, interval = 1000/30)
        anim.save(os.path.join(out_path, 'reg_anat_sample.gif'),fps = 30)
    
    return ops

def get_shift(reg_path, stack_shape, props):
    #load the ops and the registration file from the data path
    ops = np.load(os.path.join(reg_path, "ops.npy"), allow_pickle=True).flatten()[0]
    reg = np.load(os.path.join(reg_path, "reg.npy"), allow_pickle=True).flatten()[0]

    #determine the shape of the reference stack
    # stack,_ = suite2p.io.tiff.open_tiff(stack_path,sktiff=0)
    # stack_shape = props['stack_shape']

    #get the pixelwise nonrigid offsets
    xblocks, yblocks, nblocks, *_ = suite2p.registration.register.nonrigid.make_blocks(*stack_shape, ops['block_size'])
    
    xblocks = [b/stack_shape[1]*props['frame_shape'][1] for b in xblocks]
    yblocks = [b/stack_shape[0]*props['frame_shape'][0] for b in xblocks]

    x,y = suite2p.registration.register.nonrigid.upsample_block_shifts(
        props['frame_shape'][1], 
        props['frame_shape'][0], 
        nblocks, 
        xblocks, 
        yblocks, 
        reg['y_nr'] / stack_shape[0] * props['frame_shape'][0], 
        reg['x_nr'] / stack_shape[1] * props['frame_shape'][1]
    )

    #correct for the rescaling and apply rigid offsets
    x = x * stack_shape[1] / props['frame_shape'][1] - reg['x_rigid'][:,None,None]
    y = y * stack_shape[0] / props['frame_shape'][0] - reg['y_rigid'][:,None,None]
    #TODO: not sure about sign for nonrigid offset... could check phasecorr

    
    x += np.arange(0,stack_shape[1],stack_shape[1] / props['frame_shape'][1])[None,None,:]
    y += np.arange(0,stack_shape[0],stack_shape[0] / props['frame_shape'][0])[None,:,None]
    # xi = np.round(x).astype(int) + np.arange(0,stack_shape[1],stack_shape[1] / props['frame_shape'][1])[None,None,:]
    # yi = np.round(y).astype(int) + np.arange(0,stack_shape[0],stack_shape[0] / props['frame_shape'][0])[None,:,None]
    # ti = np.tile(np.arange(0,nframes)[:,None,None],(1,*props['frame_shape']))

    # mi = np.ravel_multi_index((ti,yi,xi), (nframes,*stack_shape), mode='clip')
    # O = np.zeros((nframes * stack_shape[0] * stack_shape[1]), dtype=np.int16)
    # O[mi] = movie
    # O.reshape((nframes, stack_shape[0], stack_shape[1]), inplace=True)

    #TODO: z coordinates?
    return x,y

if __name__ == "__main__":
    register(r"C:\Users\zfj\data\100522B\region1_00001.tif",r"C:\Users\zfj\data\100522B\region1_stack_00002.tif",r"C:\Users\zfj\data\100522B\region1_00001")