#copied from kb_bowtie2 on github 2019-01-16

import os.path
import subprocess


class Bowtie2Runner:

    BOWTIE2_PATH = '/kb/module/bowtie2-bin'

    def __init__(self, scratch_dir):
        self.scratch_dir = scratch_dir
        self.valid_commands = ['bowtie2',
                               'bowtie2-align-l',
                               'bowtie2-align-s',
                               'bowtie2-build',
                               'bowtie2-build-l',
                               'bowtie2-build-s',
                               'bowtie2-inspect',
                               'bowtie2-inspect-l',
                               'bowtie2-inspect-s']

    def run(self, command, options, cwd=None):
        ''' options is an array of command-line parameters passed to the RQCFilter App '''
        if command not in self.valid_commands:
            raise ValueError('Invalid bowtie2 command: ' + str(command))

        command = [os.path.join(self.BOWTIE2_PATH, command)] + options

        print('In working directory: ' + ' '.join(command))
        print('Running: ' + ' '.join(command))


        if not cwd:
          cwd = self.scratch_dir

        p = subprocess.Popen(command, cwd=cwd, shell=False)
        exitCode = p.wait()

        if (exitCode == 0):
            print('Success, exit code was: ' + str(exitCode))
        else:
            raise ValueError('Error running command: ' + ' '.join(command) + '\n' +
                             'Exit Code: ' + str(exitCode))
