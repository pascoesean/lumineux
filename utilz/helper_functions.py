from tensorly.decomposition import non_negative_parafac, non_negative_tucker
import tensorly as tl
import numpy as np


def calc_mae_for_rank(dataa, dataaa, rank, use_tucker = False):
    # calculates mean absolute error of reconstruction
    # do factorization
    # TO DO: un hard-code 'logmfi'

    nnp_weights , nnp_factors = non_negative_parafac(
        dataa, rank = rank, random_state=1
        )

    if use_tucker:
        tensor_mu, error_mu = non_negative_tucker(dataa, rank = rank, random_state=1)

        mae  = np.sum(np.abs(dataa - tl.tucker_to_tensor(tensor_mu)))
        return mae


    # get reconstruction
    reconstructed = np.sum(tl.tenalg.khatri_rao(nnp_factors, weights=nnp_weights), axis=1)

    # calc mae
    mae = np.sum(np.abs((dataaa['logmfi'].to_numpy() - reconstructed))) / len(dataaa['logmfi'].to_numpy())

    return mae
