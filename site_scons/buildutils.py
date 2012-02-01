import os
from os.path import join as pjoin
import sys
import textwrap
import re
import types
import shutil

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
    if isinstance(value, types.StringTypes):
        return value.replace(',', ' ').split()
    else:
        # Already a sequence. Return as a list
        return list(value)


def removeFile(name):
    """ Remove file (if it exists) and print a log message """
    if os.path.exists(name):
        print 'Removing file "%s"' % name
        os.remove(name)


def removeDirectory(name):
    """ Remove directory recursively and print a log message """
    if os.path.exists(name):
        print 'Removing directory "%s"' % name
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
