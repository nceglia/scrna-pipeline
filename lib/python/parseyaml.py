from genemarkermatrix import GeneMarkerMatrix
import sys

yaml = sys.argv[1]
csv  = sys.argv[2]

mat = GeneMarkerMatrix.read_yaml(yaml)
mat.write_matrix(csv)