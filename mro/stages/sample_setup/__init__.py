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
    header = samples.pop(0).split(",")
    for sample in samples:
        sample = dict(zip(header,sample.split(",")))
        outs.samples[sample["sample"]] = sample["path"]
        if sample["batch"] not in outs.batch:
            outs.batch[sample["batch"]] = []
        outs.batch[sample["batch"]].append(sample["sample"])
        if "group" in sample:
            if "group" not in outs.batch:
                outs.batch["group"] = dict()
            outs.batch["group"][sample["sample"]] = sample["group"]
            