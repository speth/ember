from __future__ import print_function
import os
from os.path import join as pjoin
import sys
import textwrap
import re
import types
import shutil
import subprocess
import platform
import SCons.Node.FS

def quoted(s):
    """ Returns the given string wrapped in double quotes."""
    return '"%s"' % s


def mglob(env, subdir, *args):
    """
    Each arg in args is assumed to be file extension,
    unless the arg starts with a '^', in which case the remainder
    of the arg is taken to be a complete pattern.
    """
    matches = []
    for ext in args:
        if ext.startswith('^'):
            matches += env.Glob(pjoin(subdir, ext[1:]))
        else:
            matches += env.Glob(pjoin(subdir, '*.%s' % ext))
    return matches


def psplit(s):
    """
    Split a path given as a string into a list.
    This is the inverse of os.path.join.
    """
    head, tail = os.path.split(s)
    path = [tail]
    while head:
        head, tail = os.path.split(head)
        path.append(tail)

    path.reverse()
    return path


def stripDrive(s):
    """
    Remove a Windows drive letter specification from a path.
    """
    if len(s) > 1 and s[1] == ':':
        return s[2:]
    else:
        return s


def which(program):
    """ Replicates the functionality of the 'which' shell command """
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

optionWrapper = textwrap.TextWrapper(initial_indent='    ',
                                   subsequent_indent='    ',
                                   width=72)

def formatOption(env, opt):
    """
    Print a nicely formatted description of a SCons configuration
    option, it's permitted values, default value, and current value
    if different from the default.
    """
    # Extract the help description from the permitted values. Original format
    # is in the format: "Help text ( value1|value2 )" or "Help text"
    if opt.help.endswith(')'):
        parts = opt.help.split('(')
        help = '('.join(parts[:-1])
        values = parts[-1][:-1].strip().replace('|', ' | ')
        if values == '':
            values = 'string'
    else:
        help = opt.help
        values = 'string'

    # Fix the representation of boolean options, which are stored as
    # Python bools, but need to be passed by the user as strings
    default = opt.default
    if default is True:
        default = 'yes'
    elif default is False:
        default = 'no'

    # First line: "* option-name: [ choice1 | choice2 ]"
    lines = ['* %s: [ %s ]' % (opt.key, values)]

    # Help text, wrapped and idented 4 spaces
    lines.extend(optionWrapper.wrap(re.sub(r'\s+', ' ',help)))

    # Default value
    lines.append('    - default: %r' % default)

    # Get the actual value in the current environment
    if opt.key in env:
        actual = env.subst('${%s}' % opt.key)
    else:
        actual = None

    # Fix the representation of boolean options
    if actual == 'True':
        actual = 'yes'
    elif actual == 'False':
        actual = 'no'

    # Print the value if it differs from the default
    if actual != default:
        lines.append('    - actual: %r' % actual)
    lines.append('')

    return lines


def listify(value):
    """
    Convert an option specified as a string to a list.  Allow both
    comma and space as delimiters. Passes lists transparently.
    """
    if isinstance(value, (list, tuple)):
        # Already a sequence. Return as a list
        return list(value)
    else:
        # assume `value` is a string
        return value.split()


def removeFile(name):
    """ Remove file (if it exists) and print a log message """
    if os.path.exists(name):
        print('Removing file "%s"' % name)
        os.remove(name)


def removeDirectory(name):
    """ Remove directory recursively and print a log message """
    if os.path.exists(name):
        print('Removing directory "%s"' % name)
        shutil.rmtree(name)


def ipdb():
    """
    Break execution and drop into an IPython debug shell at the point
    where this function is called.
    """
    from IPython.core.debugger import Pdb
    from IPython.core import ipapi

    ip = ipapi.get()
    def_colors = ip.colors
    Pdb(def_colors).set_trace(sys._getframe().f_back)


def getCommandOutput(cmd: str, *args: str, ignore_errors=False):
    """
    Run a command with arguments and return its output.
    """
    environ = dict(os.environ)
    if "PYTHONHOME" in environ:
        # Can cause problems when trying to run a different Python interpreter
        del environ["PYTHONHOME"]
    kwargs = {}
    if ignore_errors:
        kwargs["stderr"] = subprocess.DEVNULL
    data = subprocess.run(
        [cmd] + list(args),
        env=environ,
        stdout=subprocess.PIPE,
        universal_newlines=True,
        check=not ignore_errors,
        **kwargs,
    )
    return data.stdout.strip()


class DefineDict(object):
    """
    A dictionary-like object which generates appropriate preprocessor
    define statements from its dict of variable / value
    pairs. Variables whose value is None or that are not in the dict
    are left undefined.
    """
    def __init__(self, data):
        self.data = data
        self.undefined = set()

    def __getitem__(self, key):
        if key not in self.data:
            self.undefined.add(key)
            return '/* #undef %s */' % key
        elif self.data[key] is None:
            return '/* #undef %s */' % key
        else:
            return '#define %s %s' % (key, self.data[key])


class ConfigBuilder(object):
    """
    Used along with DefineDict to generate a customized config.h file
    from a config.h.in file using the variables given in 'defines'.
    """
    def __init__(self, defines):
        self.defines = DefineDict(defines)

    def __call__(self, source, target, env):
        for s, t in zip(source, target):
            config_h_in = open(str(s), "r")
            config_h = open(str(t), "w")

            config_h.write(config_h_in.read() % self.defines)
            config_h_in.close()
            config_h.close()
            self.print_config(str(t))

    def print_config(self, filename):
        print('Generating %s with the following settings:' % filename)
        for key, val in sorted(self.defines.data.items()):
            if val is not None:
                print("    %-35s %s" % (key, val))
        for key in sorted(self.defines.undefined):
            print("    %-35s %s" % (key, '*undefined*'))


# Monkey patch for SCons Cygwin bug
# See http://scons.tigris.org/issues/show_bug.cgi?id=2664
if 'cygwin' in platform.system().lower():
    SCons.Node.FS._my_normcase = lambda x: x
