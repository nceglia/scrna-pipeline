import subprocess


class Container(object):

    def __init__(self):
        self.runtime = None
        self.container_command = None

    def set_runtime(self, runtime):
        if runtime == "singularity":
            self.container_command = "singularity exec --bind /juno/work/shah/:/juno/work/shah/ --bind /work/shah/:/work/shah".split()
            self.runtime = runtime
        elif runtime == "docker":
            self.container_command = "docker run --mount type=bind,source=/Users/ceglian/,target=/Users/ceglian/ -w $(pwd)"
            self.runtime = runtime
        else:
            self.container_command = ""
            self.runtime = "local"

    def set_image(self, image):
        self.image = image

    def run(self, script, args, outs):
        cmd = [script.interpreter, script.path]
        for position, name in script.args.items():
            cmd.insert(position+2, str(getattr(args, name)))
        for position, name in script.outs.items():
            cmd.insert(position+2, str(getattr(outs, name)))
        full_cmd = self.container_command + [self.image] + cmd
        subprocess.check_output(full_cmd)
