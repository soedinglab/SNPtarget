import argparse
import pickle
import numpy as np

from ngsbiotools.parsers.rao_hic_extractor import RaoHiCParser


class HiCMatrix:

    def __init__(self, matrix, resolution, shape, matrix_type="band", annotation={}):
        self._matrix = matrix
        self._resolution = resolution
        self._shape = shape
        self._matrix_type = matrix_type
        self._annotation = annotation

    @property
    def shape(self):
        return self._shape

    @property
    def matrix(self):
        return self._matrix

    @property
    def resolution(self):
        return self._resolution

    @property
    def matrix_type(self):
        return self._matrix_type

    @property
    def annotation(self):
        return self._annotation


def main():
    parser = argparse.ArgumentParser("raohic2matrix")
    parser.add_argument("observation_file")
    parser.add_argument("--norm_file")
    parser.add_argument("resolution", type=int)
    parser.add_argument("--type", choices=["band", "full"], default="band")
    parser.add_argument("--max_distance", type=int, default=2005000)
    parser.add_argument("output_matrix")

    args = parser.parse_args()

    with open(args.observation_file) as obs_file:
        if not args.norm_file:
            hic_parser = RaoHiCParser(obs_file, args.resolution)
            ia_matrix = hic_parser.read()
        else:
            with open(args.norm_file) as norm_file:
                hic_parser = RaoHiCParser(obs_file, args.resolution, norm_file)
                ia_matrix = hic_parser.read()

    if args.type == "full":
        matrix_obj = HiCMatrix(ia_matrix, args.resolution, ia_matrix.shape,
                               matrix_type="full")
    else:
        ndiag = args.max_distance // args.resolution + 1
        diag_dict = matrix_to_diagdict(ia_matrix, ndiag)
        matrix_obj = HiCMatrix(diag_dict, args.resolution, ia_matrix.shape)

    with open(args.output_matrix, "wb") as output_file:
        pickle.dump(matrix_obj, output_file)


def matrix_to_diagdict(matrix, ndiag):
    m, n = matrix.shape
    l = min(m, n)
    ind = np.arange(l)
    diag_dict = {}
    for d in range(min(n, ndiag + 1)):
        mask = ind[0:l] + d < n
        diag_dict[d] = matrix[ind[0:l][mask], (ind[0:l]+d)[mask]].tocsr()
    return diag_dict


if __name__ == "__main__":
    main()
