import os

__MRO__ = """
stage SAMPLE_SETUP(
    in  path inventory,
    out map samples,
    out map batch,
    src py "stages/sample_setup",
)
"""

def main(args, outs):
    outs.samples = dict()
    outs.batch = dict()
    samples = open(args.inventory, "r").read().splitlines()
    samples.pop(0)
    for sample in samples:
        sample_id, batch_id, matrix = sample.split(",")
        outs.samples[sample_id] = matrix
        if batch_id not in outs.batch:
            outs.batch[batch_id] = []
        outs.batch[batch_id].append(sample_id)
            