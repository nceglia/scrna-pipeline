import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from interface.singlecellexperiment import SingleCellExperiment

class DropletUtils(object):

    def __init__(self):
        self.namespace = importr("DropletUtils")
        for attr, reference in self.namespace.__dict__.items():
            if hasattr(self,attr):
                setattr(self, "_{}".format(attr), reference)
            else:
                setattr(self, attr, reference)

    @staticmethod
    def read10xCounts(path, output):
        utils = DropletUtils()
        counts = utils.read10xCounts(path)
        sce = SingleCellExperiment.fromRS4(counts)
        sce.save(output)
        return sce



if __name__ == '__main__':
    pass
