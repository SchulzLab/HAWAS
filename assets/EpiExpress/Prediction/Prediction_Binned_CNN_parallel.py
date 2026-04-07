import json
import os
import pandas as pd
import numpy as np
import sys
import random
import re
import fnmatch
import multiprocessing

# keras is intentionally not imported here at module level.
# it must be imported inside each worker after setting thread-limiting
# env vars, otherwise TensorFlow (keras) initializes with all available cores

# Set a random seed for the results to be reproducible
seed_value = 12345
random.seed(seed_value)

def fn_patternmatch(pattern):
    """
    Grabs all files in the file system that match the pattern and returns a dictionary with {file: wildcard}.
    Only works if the wildcard is in the file name and not in a directory name.
     """
    parent_folder = '/'.join(pattern.split('/')[:-1]) + '/'
    children_pattern = pattern.split('/')[-1]
    re_pattern = re.compile(children_pattern.replace('*', '(.*)'))
    matched_files = {parent_folder + x: re_pattern.search(x).group(1)
                     for x in os.listdir(parent_folder) if fnmatch.fnmatch(x, children_pattern)}
    return matched_files


def worker_init():
    """Runs once per worker process. Sets thread limits, then imports keras."""
    # these must be set before importing keras
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["KMP_BLOCKTIME"] = "0"

    import keras  

    random.seed(seed_value)
    np.random.seed(seed_value)


def predict_gene(gene_info):
    gene, input_file, model_path, mn, mx, cnn_out = gene_info

    try:
        import keras
        model = keras.models.load_model(model_path)
        data = pd.read_csv(input_file, sep="\t", header=0)
        data = data.to_numpy()
        row_names = data[:, 0]
        x = np.log2((1 + data[:, 1:data.shape[1]].astype(float)))

        # Generate predictions and backscale the values using the given min and max.
        prediction = model.predict(x, verbose=0)
        backscaled_prediction = 2**(prediction * (mx - mn) + mn) - 1

        df_res = pd.DataFrame(prediction, index=row_names, columns=['Prediction'])
        df_res['Backscaled_prediction'] = backscaled_prediction
        output_filename = os.path.join(cnn_out, gene+'_predictions.csv.gz')
        df_res.to_csv(output_filename, sep=',', header=True, index=True, index_label='Sample')

    except Exception as e:
        print(f"Error processing gene {gene}: {e}")


def main(args):
    cores = int(args.get('cores', 1))

    output_folder = args.get("out_folder")
    model_folder = args.get("binned_cnn_model_folder")
    provided_input = args.get("provided_input")

    cnn_out = os.path.join(output_folder, 'Binned_CNN_Result')
    if not os.path.isdir(cnn_out):
        os.mkdir(cnn_out)

    input_folder = os.path.join(output_folder, 'Binned_input')
    input_files = fn_patternmatch(input_folder + '/*.txt.gz')

    model_files = {x.split('/')[-1].split('CNN_')[1].split('_w1MB')[0]: os.path.join(model_folder, x)
                   for x in os.listdir(model_folder) if x.endswith('.h5')}

    minmax_dict = pd.read_table(os.path.join(provided_input, 'MinMax_Binned_CNN.txt.gz'), sep='\t', header=0,
                                index_col=0).to_dict(orient='index')

    gene_args = []
    for input_file, gene in input_files.items():
        if gene in model_files:
            mn = minmax_dict[gene]['min']
            mx = minmax_dict[gene]['max']
            gene_args.append((gene, input_file, model_files[gene], mn, mx, cnn_out))
        else:
            print("File for gene: " + gene + " not found")

    # this ensures env vars set in worker_init take effect before keras initializes
    ctx = multiprocessing.get_context("forkserver")
    with ctx.Pool(processes=cores, initializer=worker_init) as pool:
        pool.map(predict_gene, gene_args)

    print("Binned-CNN Prediction done")


if __name__ == "__main__":
    json_file = sys.argv[1]

    try:
        with open(json_file, 'r') as file:
            args = json.load(file)
    except Exception as e:
        print(f"Error reading JSON file: {e}")
        sys.exit(1)

    main(args)