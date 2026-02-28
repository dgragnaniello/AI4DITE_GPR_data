import os, sys
from scipy.io import savemat
from paired_GPR_dataset import read_all_paired_traces, DATASET_ROOT_DIR, DATASET_SITES, MODALITIES

output_dir = "./dataset_mat/" if len(sys.argv) == 1 else sys.argv[1]
all_paired_traces = read_all_paired_traces()
for folder in DATASET_SITES:
    for i, modality in enumerate(MODALITIES):
        savemat(os.path.join(output_dir, folder, modality + '.mat'),
            {
                "original": all_paired_traces[folder][i][0],
                "aligned": all_paired_traces[folder][i][1],
            })

