import numpy as np

DATASET_ROOT_DIR = './Zenodo'
DATASET_SITES = ["Site_1", "Site_2", "Site_3"]
MODALITIES = ("ground_cart", "UAV_drone")

#EXTENSIONS_ALL = ('.shx', '.cpg', '.prj', '.kmz', '.shp', '.qmd', '.sgy', '.dbf')
EXTENSIONS_SEGY = ('.sgy', )
EXTENSIONS_TO_EXCLUDE = (".qmd", ".cpg", ".prj")


def read_segy(filename, verbose=1):  # endian '>' is Standard Rev 1 big-endian, '<' is Non-standard Rev 1 little-endian
    from segpy.reader import create_reader

    with open(filename, 'rb') as segy_in_file:
        # The seg_y_dataset is a lazy-reader, so keep the file open throughout.
        seg_y_dataset = create_reader(segy_in_file)
        if verbose>0: print(filename, 'traces:', seg_y_dataset.num_traces())

        traces = [np.asarray(seg_y_dataset.trace_samples(i)) for i in seg_y_dataset.trace_indexes()]
        traces = np.stack(traces, 1)
        try:
            headers = [np.asarray(seg_y_dataset.trace_header(i)) for i in seg_y_dataset.trace_indexes()]
        except:
            headers = None
            if verbose>0: print("Error reading trace header.")
        seg_y_dataset._fh.close()
        return traces

def read_folder(folder, modalities=MODALITIES, verbose=0):
    import os
    from geotable import load as geotable_load

    positions = {}
    traces = {}
    data_segy = {}
    data_geot = {}
    for modality in modalities:
        data_segy[modality] = {}
        data_geot[modality] = {}
        for root, dirs, files in os.walk(os.path.join(folder, modality)):
            if root.endswith("/.segpy"): continue  # cache directory

            # Loop to read files
            for file in files:
                filepath = os.path.join(root, file)

                # Skipping SEGY files, will be read later
                if file.lower().endswith(EXTENSIONS_SEGY):
                    continue

                elif not file.lower().endswith(EXTENSIONS_TO_EXCLUDE):
                    # Reading GeoTable files
                    try:
                        data_geot[modality][file] = geotable_load(filepath)
                        if verbose>1: print("read", filepath)

                    except:
                        if verbose>0: print("* error", filepath)

                else:
                    pass  # if verbose>0: print("    skipping", filepath)

        # Check that info in all GeoTable file correspond
        xy = None
        for i, geot_ in enumerate(data_geot[modality].values()):
            if verbose>2:
                for c in geot_.columns:
                    print(i, modality, c, geot_[c])
            if i == 0:
                filenames = geot_["Name"].values
                xy = [_.xy for _ in geot_['geometry_object'].values]
            else:
                assert((filenames == geot_["Name"].values).all())
                assert(xy == [_.xy for _ in geot_['geometry_object'].values])

        # Read SEGY file in the order of the GeoTable data
        max_num_samples = -1
        for filename in filenames:
            filepath = os.path.join(folder, modality, filename + ".SGY")
            data_segy[modality][filename] = read_segy(filepath)
            max_num_samples = max(max_num_samples, data_segy[modality][filename].shape[0])
            if verbose > 1: print("read", filepath)

        # Concatenate traces according to the order in the GeoTable files
        ordered_traces = []
        ordered_positions = []
        for filename, xy_ in zip(filenames, xy):
            trace = data_segy[modality][filename]
            trace = np.pad(trace, ((0, max_num_samples - len(trace)), (0, 0)), 'constant')
            ordered_traces.append(trace)

            n_ = trace.shape[1]
            xy_ = np.stack(xy_, 1)
            # Introduce a small perturbation to discriminate consecutive traces having the very same position
            same_xy = (xy_[:-1] == xy_[1:]).all(1)
            if same_xy.any():
                delta = np.diff(xy_, axis=0).mean(0) / 10.
                inds = same_xy.nonzero()[0]
                for ind in inds:
                    xy_[ind+1] = xy_[ind] + delta
            xy_ = np.stack([
                np.interp(np.linspace(0, xy_.shape[0] - 1, n_), np.linspace(0, xy_.shape[0] - 1, xy_.shape[0]), xy_[:, 0]),
                np.interp(np.linspace(0, xy_.shape[0] - 1, n_), np.linspace(0, xy_.shape[0] - 1, xy_.shape[0]), xy_[:, 1])
            ],1)
            ordered_positions.append(xy_)


        traces[modality] = ordered_traces
        positions[modality] = ordered_positions

        traces[modality] = np.concatenate(traces[modality],1)
        positions[modality] = np.concatenate(positions[modality],0)

    return positions, traces, data_segy, data_geot

def compute_modalities_correspondance(position_modality_A, position_modality_B):
    from scipy.spatial.distance import cdist

    D = cdist(position_modality_A, position_modality_B)
    modB2modA = D.argmin(axis=1)
    modA2modB = D.argmin(axis=0)
    return modB2modA, modA2modB

def get_paired_traces(traces_modality_A, traces_modality_B, modB2modA, index=None, plot=False):
    from optype.inspect import is_iterable

    lA, nA = traces_modality_A.shape
    lB, nB = traces_modality_B.shape
    if index is None:
        index = range(nA)
    elif not is_iterable(index):
        index = [index]

    traces_B = np.zeros((lB, nA))
    for i in index:
        traces_B[:, i] = traces_modality_B[:, modB2modA[i]]

    if plot:
        import matplotlib.pyplot as plt
        tA = np.arange(lA)
        tB = np.arange(lB)
        for i in index:
            plt.figure()
            plt.plot(tA, traces_modality_A[:, i])
            plt.plot(tB, traces_B[:, i])
            plt.title("Mod. A no. %d --> Mod B no. %d" % (i, modB2modA[i]))
            plt.legend(["A", "B"])
            plt.show()
    return traces_B

def plot_geot_layers(layers, names=None, title=None):
    import matplotlib.pyplot as plt
    from shapely import plotting as shpl

    for layer in layers:
        shpl.plot_points(layer.draw(), markersize=2)

    if names is None:
        names = ["layer #%d" % i for i, layer in enumerate(layers)]
    plt.legend(names)
    if title is not None:
        plt.title(title)

def read_paired_traces(folder, verbose=1, plot=False, modality_A=MODALITIES[0], modality_B=MODALITIES[1]):
    positions, traces, data_segy, data_geot = read_folder(folder, verbose=verbose)
    modB2modA, modA2modB = compute_modalities_correspondance(positions[modality_A], positions[modality_B])
    traces_modB2modA = get_paired_traces(traces[modality_A], traces[modality_B], modB2modA, plot=plot)
    traces_modA2modB = get_paired_traces(traces[modality_B], traces[modality_A], modA2modB, plot=plot)
    return (traces[modality_A], traces_modB2modA), (traces[modality_B], traces_modA2modB)

def read_all_paired_traces(verbose=1, plot=False):
    pared_traces = {}
    for folder in DATASET_SITES:
        folder = os.path.join(DATASET_ROOT_DIR, folder)
        pared_traces[folder] = read_paired_traces(folder, verbose, plot)
    return pared_traces



if __name__ == "__main__":
    import os, sys
    plot = True if len(sys.argv) == 1 else int(sys.argv[1]) > 0

    positions = {}
    traces = {}
    data_segy = {}
    data_geot = {}
    for folder in DATASET_SITES:
        folder = os.path.join(DATASET_ROOT_DIR, folder)
        # Eventually extract the zip file
        if not os.path.isdir(folder):
            import zipfile
            with zipfile.ZipFile(folder + ".zip", 'r') as zip_ref:
                zip_ref.extractall(DATASET_ROOT_DIR)

        # Read the site data
        positions[folder], traces[folder], data_segy[folder], data_geot[folder] = read_folder(folder, verbose=1)

        if plot:
            import matplotlib.pyplot as plt

            modB2modA, modA2modB = compute_modalities_correspondance(positions[folder][MODALITIES[0]], positions[folder][MODALITIES[1]])
            traces_modB2modA = get_paired_traces(traces[folder][MODALITIES[0]], traces[folder][MODALITIES[1]], modB2modA, plot=False)
            traces_modA2modB = get_paired_traces(traces[folder][MODALITIES[1]], traces[folder][MODALITIES[0]], modA2modB, plot=False)
            v_ = 15000


            plt.figure(figsize=(20,15))
            plt.subplot(1,2,1)
            plot_geot_layers(list(data_geot[folder][MODALITIES[0]].values()), title="Folder: {}, Modality: {}".format(os.path.basename(folder), MODALITIES[0]))
            plt.subplot(1,2,2)
            plot_geot_layers(list(data_geot[folder][MODALITIES[1]].values()), title="Folder: {}, Modality: {}".format(os.path.basename(folder), MODALITIES[1]))

            plt.figure(figsize=(20,15))
            plt.subplot(4,1,1)
            plt.imshow(traces[folder][MODALITIES[0]], vmin=-v_, vmax=v_)
            plt.title("Folder: {}, Modality: {}".format(os.path.basename(folder), MODALITIES[0]))

            plt.subplot(4,1,2)
            plt.imshow(traces[folder][MODALITIES[1]], vmin=-v_, vmax=v_)
            plt.title(MODALITIES[1])

            plt.subplot(4,1,3)
            plt.imshow(traces_modB2modA, vmin=-v_, vmax=v_)
            plt.title("Folder: {}, from {} to {}".format(os.path.basename(folder), MODALITIES[1],MODALITIES[0]))

            plt.subplot(4,1,4)
            plt.imshow(traces_modA2modB, vmin=-v_, vmax=v_)
            plt.title("Folder: {}, from {} to {}".format(os.path.basename(folder), MODALITIES[0],MODALITIES[1]))


            plt.figure(figsize=(20,15))
            trace_index_to_plot = 100
            plt.subplot(2,1,1)
            plt.title("Folder: {}, Modality: {} - Traces n. {:d}".format(os.path.basename(folder), MODALITIES[0], trace_index_to_plot))
            plt.plot(traces[folder][MODALITIES[0]][:, trace_index_to_plot])
            plt.plot(traces_modB2modA[:, trace_index_to_plot])
            plt.legend((MODALITIES[0], "from {} to {}".format(MODALITIES[1],MODALITIES[0])))

            trace_index_to_plot = 100
            plt.subplot(2,1,2)
            plt.title("Folder: {}, Modality: {} - Traces n. {:d}".format(os.path.basename(folder), MODALITIES[1], trace_index_to_plot))
            plt.plot(traces_modA2modB[:, trace_index_to_plot])
            plt.plot(traces[folder][MODALITIES[1]][:, trace_index_to_plot])
            plt.legend(("from {} to {}".format(MODALITIES[0], MODALITIES[1]), MODALITIES[1]))

            plt.show()

    print("Done.")
