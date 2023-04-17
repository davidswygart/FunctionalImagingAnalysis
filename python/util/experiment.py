if __name__ == "__main__":
    import preprocess
else:
    from . import preprocess

import numpy as np


#TODO: we should accept a list of raw files
#we can save the flips as an offset from the first file, maybe?
#instead of saving the frame index for every pixel, we can just save the index for each flip?

#TODO: separate out the common experiment code from the chirp/spot/bar code?

class Experiment:
    def __init__(self, timing_file: str, raw_file: str, bin_path: str):
        self.timing_file = timing_file
        self.raw_file = raw_file
        self.bin_path = bin_path

        self.epochs = preprocess.flatten_epochs(preprocess.extract_parameters(self.timing_file))
        self.epochs['LED'] = (self.epochs.trial_type == b'field') * self.epochs.spotLED + (self.epochs.trial_type == b'chirp') * self.epochs.chirpLED + (self.epochs.trial_type == b'bars') * self.epochs.barLED
        self.epochs['end_time'] = np.uint64(self.epochs.preTime + self.epochs.stimTime + self.epochs.tailTime) * 10000 + self.epochs.start_time
        self.epochs = self.epochs.sort_values('start_time')

        #TODO: we should accept a list of raw files... if props are the same, then combine them
        # in doing so, we need to make sure we're registering to the correct region...
        self.raw, self.props = preprocess.read_raw_file(self.raw_file, self.bin_path)
    

        # the trial times are encoded in the epochs as well as the trigger channel
        # let's run edge detection on the trigger channel first
        
        self.trigger, self.flips = preprocess.get_trigger_times(self.raw, self.props)

        # we can then segment the epochs and cluster the trigger channel into 4 groups
        # NOTE: group the epochs by the LED/NDF settings, since this affects the trigger level

        #{pseudocode}
        # epochs.groupby(['LED','NDF']).apply( lambda epoch:
        #    first_ts = np.argwhere(props['time_stamps']<epoch.start_time)[-1]
        #    last_ts = np.argwhere(props['time_stamps']>epoch.end_time)[0]
        #    these_flips.append...
        # )


#TODO: triggering is done by flipping every other line so that time points are sequential
# we have to be careful when using this, since the trigger times don't index into the raw image
# we may want to flip them back before using them, but not sure when to do this? it makes it harder
# to figure out exactly which trigger corresponds to a pixel