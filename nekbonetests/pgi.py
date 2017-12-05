#!/usr/bin/env python2

import os, re, stat, sys, unittest
from subprocess import call, PIPE, Popen, STDOUT

class NekboneTestCase(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        self.source_root = os.path.join(
            os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
            'src'
        )

        unittest.TestCase.__init__(self, *args, **kwargs)

    def config_makenek(self, opts, infile, outfile):
        """ Redefine variables in a makenek file

        Given a path to infile, redefine the variables & values in opts, then
        output to outfile.  The infile and outfile can be the same file.

        Only the options already available in makenek may be set.  This function
        doesn't allow you to declare arbitrary variables in makenek.
        (TODO: Raise warning when attempting to set an unavailable variable.)

        Args:
            opts ({variable : value}): Set each "variable=value" in makenek
            infile (): The input makenek file
            outfile (): The output makenek file

        """
        with open(infile, 'r') as f:
            lines = f.readlines()

        for key, val in opts.iteritems():
            lines = [re.sub(r'^#*{0}=\"+.*?\"+'.format(key), r'{0}="{1}"'.format(key, val), l) for l in lines]

        lines = [re.sub(r'(^source\s+\$SOURCE_ROOT/makenek.inc)', r'\g<1> >compiler.out', l)
                 for l in lines]

        lines = [re.sub(r'(.+)2>&1\s+\|\s*tee\s+compiler.out', r'\g<1>', l)
                 for l in lines]

        with open(outfile, 'w') as f:
            f.writelines(lines)
        os.chmod(outfile,
                 stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH |
                 stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH |
                 stat.S_IWUSR)

    def build_nekbone(self, source_root, test_root, opts=None, verbose=False):

        if not opts:
            _opts = {}
        else:
            _opts = opts.copy()
        _opts.update(SOURCE_ROOT=source_root)

        print('Compiling nek5000...')
        print('    Using source directory "{0}"'.format(source_root))
        print('    Using working directory "{0}"'.format(test_root))
        for key, val in _opts.iteritems():
            print('    Using {0}="{1}"'.format(key, val))

        makenek     = os.path.join(test_root, 'makenek')
        logfile     = os.path.join(test_root, 'compiler.out')
        try:
            self.config_makenek(
                opts=_opts,
                infile=makenek,
                outfile=makenek
            )

            call([makenek, 'clean'], cwd=test_root)
            if verbose:
                with open(logfile, 'w') as f:
                    proc = Popen([makenek], cwd=test_root, stderr=STDOUT, stdout=PIPE)
                    for line in proc.stdout:
                        sys.stdout.write(line)
                        f.write(line)
            else:
                with open(logfile, 'w') as f:
                    call([makenek], cwd=test_root, stdout=f)

        except:
            print('Could not compile nek5000!')
            raise
        else:
            print('Successfully compiled nek5000!')

    @staticmethod
    def run_nek(cwd, rea_file, ifmpi, log_suffix='', n_procs=1, step_limit=None, verbose=False):
        # Paths to executables, files
        nek5000      = os.path.join(cwd, 'nek5000')
        logfile      = os.path.join(cwd, '{0}.log.{1}{2}'.format(rea_file, n_procs, log_suffix))
        session_name = os.path.join(cwd, 'SESSION.NAME')
        ioinfo       = os.path.join(cwd, 'ioinfo')
        sch_file     = os.path.join(cwd, '{0}.sch'.format(rea_file))
        if ifmpi:
            command = ['mpiexec', '-np', str(n_procs), nek5000]
        else:
            command = [nek5000]

        print("Running nek5000...")
        print('    Using command "{0}"'.format(' '.join(command)))
        print('    Using working directory "{0}"'.format(cwd))
        print('    Using .rea file "{0}"'.format(rea_file))

        # An OSError here can be expected
        # If the examples directory is clean, there will be no .sch file and
        # os.remove(sch_file) will be expected to fail.
        try:
            os.remove(sch_file)
        except OSError as E:
            # TODO: Change to warnings.warning
            print("    Could not remove {0}: {1}".format(sch_file, E))

        # Any error here is unexepected
        try:
            with open(session_name, 'w') as f:
                f.writelines([rea_file+"\n", cwd+'/\n'])

            if step_limit:
                with open(ioinfo, 'w') as f:
                    f.writelines(['-{0}'.format(step_limit)])

            if verbose:
                with open(logfile, 'w') as f:
                    proc =Popen(command, cwd=cwd, stderr=STDOUT, stdout=PIPE)
                    for line in proc.stdout:
                        sys.stdout.write(line)
                        f.write(line)
            else:
                with open(logfile, 'w') as f:
                    call(command, cwd=cwd, stdout=f)

        except Exception as E:
            # TODO: Change to warnings.warn()
            print('Could not successfully run nek5000! Caught error: {0}'.format(E))
        else:
            print('Finished running nek5000!')


class Example1(NekboneTestCase):

#     def setUp(self):
#         self.config_makenek()

    def __init__(self, *args, **kwargs):
        NekboneTestCase.__init__(self, *args, **kwargs)
        self.test_root = os.path.join(
            os.path.dirname(self.source_root),
            'test',
            'example1'
        )

    def setUp(self):
        opts = dict(
            CC = "gcc",
            F77 = "gfortran",
            IFMPI = "false",
        )
        self.build_nekbone(
            source_root = self.source_root,
            test_root = self.test_root,
            opts = opts,
            verbose = True
        )

    def test(self):

        print self.source_root
        print self.test_root








