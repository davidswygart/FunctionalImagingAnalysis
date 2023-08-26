from util.experiment import Experiment

if __name__ == "__main__":
    symphony_file = r"C:\Users\zfj\data\021923B\021923B.h5"
    raw_file_root = r"p1141_p1632_00002"
    raw_file = f"C:\\Users\\zfj\\data\\021923B\\{raw_file_root}.tif"
    artifact_file = f"C:\\Users\\zfj\\code\\functionalimaginganalysis\\python\\artifact_lut.npy"
    
    e = Experiment(symphony_file, [raw_file], None, artifact_file=artifact_file)